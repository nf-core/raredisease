//
// A quality check subworkflow for processed bams.
//

include { CHROMOGRAPH as CHROMOGRAPH_COV                          } from '../../../modules/nf-core/chromograph/main'
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
include { VERIFYBAMID_VERIFYBAMID2                                } from '../../../modules/nf-core/verifybamid/verifybamid2/main'

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
        ch_svd_bed                      // channel: [optional] [ path(bed) ]
        ch_svd_mu                       // channel: [optional] [ path(meanpath) ]
        ch_svd_ud                       // channel: [optional] [ path(ud) ]
        ch_sambamba_bed                 // channel: [optional] [ val(meta), path(bed) ]
        ch_target_intervals             // channel: [mandatory] [ path(intervals_list) ]
        val_analysis_type               // string: "wes", "wgs", or "mito"
        val_aligner                     // string: "bwa", "bwamem2", "bwameme", or "sentieon"
        val_target_bed                  // string: path to target bed file
        skip_ngsbits                    // boolean

    main:
        ch_cov       = channel.empty()
        ch_cov_y     = channel.empty()
        ch_hsmetrics = channel.empty()
        ch_ngsbits   = channel.empty()

        PICARD_COLLECTMULTIPLEMETRICS (ch_bam_bai, ch_genome_fasta, ch_genome_fai)

        ch_bam_bai
            .combine(ch_bait_intervals)
            .combine(ch_target_intervals)
            .set { ch_hsmetrics_in}

        if (val_target_bed) {
            ch_hsmetrics = PICARD_COLLECTHSMETRICS (ch_hsmetrics_in, ch_genome_fasta, ch_genome_fai, [[:],[]], [[:],[]]).metrics
        }

        TIDDIT_COV (ch_bam_bai, [[],[]]) // 2nd pos. arg is req. only for cram input

        UCSC_WIGTOBIGWIG (TIDDIT_COV.out.wig, ch_genome_chrsizes)

        CHROMOGRAPH_COV([[:],[]], TIDDIT_COV.out.wig, [[:],[]], [[:],[]], [[:],[]], [[:],[]], [[:],[]])

        ch_bam_bai.map{ meta, bam, bai -> [meta, bam, bai, []]}.set{ch_mosdepth_in}
        MOSDEPTH (ch_mosdepth_in, ch_genome_fasta)


        SAMBAMBA_DEPTH(ch_bam_bai, ch_sambamba_bed, 'region')

        // COLLECT WGS METRICS
        if (!val_analysis_type.equals("wes")) {
            if (val_aligner.matches("bwa|bwameme|bwamem2")) {
                PICARD_COLLECTWGSMETRICS_WG ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_intervals_wgs )
                PICARD_COLLECTWGSMETRICS_Y ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_intervals_y )
                ch_cov      = PICARD_COLLECTWGSMETRICS_WG.out.metrics
                ch_cov_y    = PICARD_COLLECTWGSMETRICS_Y.out.metrics
            } else if (val_aligner.equals("sentieon")) {
                SENTIEON_WGSMETRICS_WG ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_intervals_wgs.map{ interval -> [[:], interval]} )
                SENTIEON_WGSMETRICS_Y ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_intervals_y.map{ interval -> [[:], interval]} )
                ch_cov      = SENTIEON_WGSMETRICS_WG.out.wgs_metrics
                ch_cov_y    = SENTIEON_WGSMETRICS_Y.out.wgs_metrics
            }
        }
        // Check sex
        if (!skip_ngsbits) {
            NGSBITS_SAMPLEGENDER(ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_ngsbits_method)
            ch_ngsbits  = NGSBITS_SAMPLEGENDER.out.tsv
        }
        // Check contamination
        ch_svd_in = ch_svd_ud.combine(ch_svd_mu).combine(ch_svd_bed).collect()
        VERIFYBAMID_VERIFYBAMID2(ch_bam_bai, ch_svd_in, [], ch_genome_fasta.map {_meta, fasta-> fasta})

        ch_qc_bam_files = PICARD_COLLECTMULTIPLEMETRICS.out.metrics.transpose()
            .mix(PICARD_COLLECTMULTIPLEMETRICS.out.pdf.transpose())
            .mix(TIDDIT_COV.out.wig)
            .mix(TIDDIT_COV.out.cov)
            .mix(UCSC_WIGTOBIGWIG.out.bw)
            .mix(CHROMOGRAPH_COV.out.plots.transpose())
            .mix(MOSDEPTH.out.global_txt)
            .mix(MOSDEPTH.out.summary_txt)
            .mix(MOSDEPTH.out.per_base_d4)
            .mix(MOSDEPTH.out.regions_txt)
            .mix(MOSDEPTH.out.per_base_bed)
            .mix(MOSDEPTH.out.per_base_csi)
            .mix(MOSDEPTH.out.regions_bed)
            .mix(MOSDEPTH.out.regions_csi)
            .mix(MOSDEPTH.out.quantized_bed)
            .mix(MOSDEPTH.out.quantized_csi)
            .mix(MOSDEPTH.out.thresholds_bed)
            .mix(MOSDEPTH.out.thresholds_csi)
            .mix(SAMBAMBA_DEPTH.out.bed)
            .mix(VERIFYBAMID_VERIFYBAMID2.out.self_sm)
            .mix(VERIFYBAMID_VERIFYBAMID2.out.log)
            .mix(VERIFYBAMID_VERIFYBAMID2.out.ud)
            .mix(VERIFYBAMID_VERIFYBAMID2.out.bed)
            .mix(VERIFYBAMID_VERIFYBAMID2.out.mu)
            .mix(VERIFYBAMID_VERIFYBAMID2.out.ancestry)
            .mix(ch_hsmetrics)
            .mix(ch_cov)
            .mix(ch_cov_y)

    emit:
        bigwig           = UCSC_WIGTOBIGWIG.out.bw                   // channel: [ val(meta), path(bw) ]
        cov              = ch_cov                                    // channel: [ val(meta), path(metrics) ]
        cov_y            = ch_cov_y                                  // channel: [ val(meta), path(metrics) ]
        d4               = MOSDEPTH.out.per_base_d4                  // channel: [ val(meta), path(d4) ]
        global_dist      = MOSDEPTH.out.global_txt                   // channel: [ val(meta), path(txt) ]
        hs_metrics       = ch_hsmetrics                              // channel: [ val(meta), path(metrics) ]
        multiple_metrics = PICARD_COLLECTMULTIPLEMETRICS.out.metrics // channel: [ val(meta), path(metrics) ]
        qc_bam_files     = ch_qc_bam_files                           // channel: [ val(meta), path(file) ]
        self_sm          = VERIFYBAMID_VERIFYBAMID2.out.self_sm      // channel: [ val(meta), path(selfSM) ]
        sex_check        = ch_ngsbits                                // channel: [ val(meta), path(tsv) ]
        tiddit_wig       = TIDDIT_COV.out.wig                        // channel: [ val(meta), path(wig) ]
}
