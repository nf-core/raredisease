//
// A quality check subworkflow for processed bams.
//

include { CHROMOGRAPH as CHROMOGRAPH_COV                          } from '../../../modules/nf-core/chromograph/main'
include { RIKER_MULTI                                            } from '../../../modules/nf-core/riker/multi/main'
include { MOSDEPTH                                                } from '../../../modules/nf-core/mosdepth/main'
include { NGSBITS_SAMPLEGENDER                                    } from '../../../modules/nf-core/ngsbits/samplegender/main'
include { PICARD_COLLECTHSMETRICS                                 } from '../../../modules/nf-core/picard/collecthsmetrics/main'
include { PICARD_COLLECTMULTIPLEMETRICS                           } from '../../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTWGSMETRICS as PICARD_COLLECTWGSMETRICS_WG } from '../../../modules/nf-core/picard/collectwgsmetrics/main'
include { PICARD_COLLECTWGSMETRICS as PICARD_COLLECTWGSMETRICS_Y  } from '../../../modules/nf-core/picard/collectwgsmetrics/main'
include { SAMBAMBA_DEPTH                                          } from '../../../modules/nf-core/sambamba/depth/main'
include { SENTIEON_WGSMETRICS as SENTIEON_WGSMETRICS_WG           } from '../../../modules/nf-core/sentieon/wgsmetrics/main'
include { SENTIEON_WGSMETRICS as SENTIEON_WGSMETRICS_Y            } from '../../../modules/nf-core/sentieon/wgsmetrics/main'
include { TIDDIT_COV                                              } from '../../../modules/nf-core/tiddit/cov/main'
include { UCSC_WIGTOBIGWIG                                        } from '../../../modules/nf-core/ucsc/wigtobigwig/main'
workflow QC_BAM {

    take:
        ch_bam_bai                      // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_bait_intervals               // channel: [mandatory] [ path(intervals_list) ]
        ch_genome_chrsizes              // channel: [mandatory] [ path(sizes) ]
        ch_genome_fai                   // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_fasta                 // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_intervals_wgs                // channel: [mandatory] [ path(intervals) ]
        ch_intervals_y                  // channel: [mandatory] [ path(intervals) ]
        ch_ngsbits_method               // channel: [val(method)]
        ch_sambamba_bed                 // channel: [optional] [ val(meta), path(bed) ]
        ch_target_intervals             // channel: [mandatory] [ path(intervals_list) ]
        val_analysis_type               // string: "wes", "wgs", or "mito"
        val_aligner                     // string: "bwa", "bwamem2", "bwameme", or "sentieon"
        val_target_bed                  // string: path to target bed file
        skip_ngsbits                    // boolean
        val_qc_metrics_tool             // string: "picard" or "riker"

    main:
        ch_cov                    = channel.empty()
        ch_cov_y                  = channel.empty()
        ch_hsmetrics              = channel.empty()
        ch_ngsbits                = channel.empty()
        ch_picard_multimetrics    = channel.empty()
        ch_picard_multimetrics_pdf = channel.empty()
        ch_riker_alignment        = channel.empty()
        ch_riker_wgs              = channel.empty()
        ch_riker_isize            = channel.empty()
        ch_riker_base             = channel.empty()
        ch_riker_mean_qual        = channel.empty()
        ch_riker_qual_dist        = channel.empty()
        ch_riker_hybcap           = channel.empty()
        ch_riker_gcbias           = channel.empty()

        // COLLECT WGS METRICS
        // The WGS coverage metrics source depends only on the aligner, not on qc_metrics_tool:
        // sentieon -> Sentieon WgsMetricsAlgo; bwa-family + picard -> Picard CollectWgsMetrics;
        // bwa-family + riker -> RIKER_MULTI (the wgs tool, emitted as riker_wgs_metrics below).
        if (val_analysis_type != 'wes') {
            if (val_aligner == 'sentieon') {
                SENTIEON_WGSMETRICS_WG ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_intervals_wgs.map{ interval -> [[:], interval]} )
                SENTIEON_WGSMETRICS_Y ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_intervals_y.map{ interval -> [[:], interval]} )
                ch_cov      = SENTIEON_WGSMETRICS_WG.out.wgs_metrics
                ch_cov_y    = SENTIEON_WGSMETRICS_Y.out.wgs_metrics
            } else if (val_qc_metrics_tool == 'picard' && val_aligner in ['bwa', 'bwameme', 'bwamem2']) {
                PICARD_COLLECTWGSMETRICS_WG ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_intervals_wgs )
                PICARD_COLLECTWGSMETRICS_Y ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_intervals_y )
                ch_cov      = PICARD_COLLECTWGSMETRICS_WG.out.metrics
                ch_cov_y    = PICARD_COLLECTWGSMETRICS_Y.out.metrics
            }
        }

        if (val_qc_metrics_tool == 'picard') {
            PICARD_COLLECTMULTIPLEMETRICS (ch_bam_bai, ch_genome_fasta, ch_genome_fai)
            ch_picard_multimetrics     = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
            ch_picard_multimetrics_pdf = PICARD_COLLECTMULTIPLEMETRICS.out.pdf

            ch_bam_bai
                .combine(ch_bait_intervals)
                .combine(ch_target_intervals)
                .set { ch_hsmetrics_in}

            if (val_target_bed) {
                ch_hsmetrics = PICARD_COLLECTHSMETRICS (ch_hsmetrics_in, ch_genome_fasta, ch_genome_fai, [[:],[]], [[:],[]]).metrics
            }
        } else if (val_qc_metrics_tool == 'riker') {
            // Single source of truth for riker tool selection; conf/modules/qc_bam.config only
            // formats meta.riker_tools into the --tools argument. riker collects 'wgs' itself only
            // on non-wes runs with a bwa-family aligner; for sentieon the WGS coverage comes from
            // WgsMetricsAlgo above, so 'wgs' is omitted here to avoid computing it twice. 'hybcap'
            // requires a target bed.
            def riker_tools = ['alignment', 'basic', 'isize', 'gcbias']
            if (val_analysis_type != 'wes' && val_aligner != 'sentieon') { riker_tools << 'wgs' }
            if (val_target_bed) { riker_tools << 'hybcap' }

            // When no target bed is supplied, ch_bait_intervals/ch_target_intervals are empty
            // channels; combining against them would yield an empty channel and RIKER_MULTI would
            // never run. Build the input with empty bait/target placeholders in that case instead.
            ch_riker_in = ( val_target_bed
                ? ch_bam_bai
                    .combine(ch_bait_intervals)
                    .combine(ch_target_intervals)
                : ch_bam_bai.map { meta, bam, bai -> [meta, bam, bai, [], []] }
            ).map { meta, bam, bai, baits, targets -> [meta + [riker_tools: riker_tools], bam, bai, baits, targets] }

            RIKER_MULTI (
                ch_riker_in,
                // .first() keeps the reference a value channel so it broadcasts to every sample;
                // a bare .join() is a queue channel and would pair with only the first sample.
                ch_genome_fasta.join(ch_genome_fai, failOnMismatch:true, failOnDuplicate:true).first()
            )

            // RIKER_MULTI also emits error_*, gcbias_detail, hybcap_per_target/per_base and pdf
            // channels; those are intentionally not wired here as the pipeline does not consume them.
            ch_riker_alignment = RIKER_MULTI.out.alignment_metrics
            ch_riker_wgs       = RIKER_MULTI.out.wgs_metrics
            ch_riker_isize     = RIKER_MULTI.out.isize_metrics
            ch_riker_base      = RIKER_MULTI.out.base_dist
            ch_riker_mean_qual = RIKER_MULTI.out.mean_qual
            ch_riker_qual_dist = RIKER_MULTI.out.qual_dist
            ch_riker_hybcap    = RIKER_MULTI.out.hybcap_metrics
            ch_riker_gcbias    = RIKER_MULTI.out.gcbias_summary
        } else {
            error "Unknown qc_metrics_tool '${val_qc_metrics_tool}'; expected 'picard' or 'riker'."
        }

        TIDDIT_COV (ch_bam_bai, [[],[]]) // 2nd pos. arg is req. only for cram input

        UCSC_WIGTOBIGWIG (TIDDIT_COV.out.wig, ch_genome_chrsizes)

        CHROMOGRAPH_COV([[:],[]], TIDDIT_COV.out.wig, [[:],[]], [[:],[]], [[:],[]], [[:],[]], [[:],[]])

        ch_bam_bai.map{ meta, bam, bai -> [meta, bam, bai, []]}.set{ch_mosdepth_in}
        MOSDEPTH (ch_mosdepth_in, ch_genome_fasta)


        SAMBAMBA_DEPTH(ch_bam_bai, ch_sambamba_bed, 'region')
        // Check sex
        if (!skip_ngsbits) {
            NGSBITS_SAMPLEGENDER(ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_ngsbits_method)
            ch_ngsbits  = NGSBITS_SAMPLEGENDER.out.tsv
        }
    emit:
        chromograph_cov_plots                 = CHROMOGRAPH_COV.out.plots                    // channel: [ val(meta), path(png) ]
        mosdepth_global_txt                   = MOSDEPTH.out.global_txt                      // channel: [ val(meta), path(txt) ]
        mosdepth_per_base_bed                 = MOSDEPTH.out.per_base_bed                    // channel: [ val(meta), path(bed.gz) ]
        mosdepth_per_base_csi                 = MOSDEPTH.out.per_base_csi                    // channel: [ val(meta), path(csi) ]
        mosdepth_per_base_d4                  = MOSDEPTH.out.per_base_d4                     // channel: [ val(meta), path(d4) ]
        mosdepth_quantized_bed                = MOSDEPTH.out.quantized_bed                   // channel: [ val(meta), path(bed.gz) ]
        mosdepth_quantized_csi                = MOSDEPTH.out.quantized_csi                   // channel: [ val(meta), path(csi) ]
        mosdepth_regions_bed                  = MOSDEPTH.out.regions_bed                     // channel: [ val(meta), path(bed.gz) ]
        mosdepth_regions_csi                  = MOSDEPTH.out.regions_csi                     // channel: [ val(meta), path(csi) ]
        mosdepth_regions_txt                  = MOSDEPTH.out.regions_txt                     // channel: [ val(meta), path(txt) ]
        mosdepth_summary_txt                  = MOSDEPTH.out.summary_txt                     // channel: [ val(meta), path(txt) ]
        mosdepth_thresholds_bed               = MOSDEPTH.out.thresholds_bed                  // channel: [ val(meta), path(bed.gz) ]
        mosdepth_thresholds_csi               = MOSDEPTH.out.thresholds_csi                  // channel: [ val(meta), path(csi) ]
        ngsbits_samplegender_tsv              = ch_ngsbits                                   // channel: [ val(meta), path(tsv) ]
        picard_collecthsmetrics_metrics       = ch_hsmetrics                                 // channel: [ val(meta), path(metrics) ]
        picard_collectmultiplemetrics_metrics = ch_picard_multimetrics                        // channel: [ val(meta), path(metrics) ]
        picard_collectmultiplemetrics_pdf     = ch_picard_multimetrics_pdf                   // channel: [ val(meta), path(pdf) ]
        sambamba_depth_bed                    = SAMBAMBA_DEPTH.out.bed                       // channel: [ val(meta), path(bed) ]
        tiddit_cov_cov                        = TIDDIT_COV.out.cov                           // channel: [ val(meta), path(bed) ]
        tiddit_cov_wig                        = TIDDIT_COV.out.wig                           // channel: [ val(meta), path(wig) ]
        ucsc_wigtobigwig_bw                   = UCSC_WIGTOBIGWIG.out.bw                      // channel: [ val(meta), path(bw) ]
        wgsmetrics_wg                         = ch_cov                                       // channel: [ val(meta), path(metrics) ]
        wgsmetrics_y                          = ch_cov_y                                     // channel: [ val(meta), path(metrics) ]
        riker_alignment_metrics               = ch_riker_alignment                           // channel: [ val(meta), path(txt) ]
        riker_wgs_metrics                     = ch_riker_wgs                                 // channel: [ val(meta), path(txt) ]
        riker_isize_metrics                   = ch_riker_isize                               // channel: [ val(meta), path(txt) ]
        riker_base_dist                       = ch_riker_base                                // channel: [ val(meta), path(txt) ]
        riker_mean_qual                       = ch_riker_mean_qual                           // channel: [ val(meta), path(txt) ]
        riker_qual_dist                       = ch_riker_qual_dist                           // channel: [ val(meta), path(txt) ]
        riker_hybcap_metrics                  = ch_riker_hybcap                              // channel: [ val(meta), path(txt) ]
        riker_gcbias_summary                  = ch_riker_gcbias                              // channel: [ val(meta), path(txt) ]
}
