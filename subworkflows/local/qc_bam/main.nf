//
// A quality check subworkflow for processed bams.
//

include { PICARD_COLLECTMULTIPLEMETRICS                            } from '../../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTHSMETRICS                                  } from '../../../modules/nf-core/picard/collecthsmetrics/main'
include { CHROMOGRAPH as CHROMOGRAPH_COV                           } from '../../../modules/nf-core/chromograph/main'
include { QUALIMAP_BAMQC                                           } from '../../../modules/nf-core/qualimap/bamqc/main'
include { TIDDIT_COV                                               } from '../../../modules/nf-core/tiddit/cov/main'
include { MOSDEPTH                                                 } from '../../../modules/nf-core/mosdepth/main'
include { SAMBAMBA_DEPTH                                           } from '../../../modules/nf-core/sambamba/depth/main'
include { UCSC_WIGTOBIGWIG                                         } from '../../../modules/nf-core/ucsc/wigtobigwig/main'
include { PICARD_COLLECTWGSMETRICS as PICARD_COLLECTWGSMETRICS_WG  } from '../../../modules/nf-core/picard/collectwgsmetrics/main'
include { PICARD_COLLECTWGSMETRICS as PICARD_COLLECTWGSMETRICS_Y   } from '../../../modules/nf-core/picard/collectwgsmetrics/main'
include { SENTIEON_WGSMETRICS as SENTIEON_WGSMETRICS_WG            } from '../../../modules/nf-core/sentieon/wgsmetrics/main'
include { SENTIEON_WGSMETRICS as SENTIEON_WGSMETRICS_Y             } from '../../../modules/nf-core/sentieon/wgsmetrics/main'
include { NGSBITS_SAMPLEGENDER                                     } from '../../../modules/nf-core/ngsbits/samplegender/main'
include { VERIFYBAMID_VERIFYBAMID2                                 } from '../../../modules/nf-core/verifybamid/verifybamid2/main'

workflow QC_BAM {

    take:
        ch_bam                          // channel: [mandatory] [ val(meta), path(bam) ]
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
        skip_qualimap                   // boolean

    main:
        ch_cov       = channel.empty()
        ch_cov_y     = channel.empty()
        ch_versions  = channel.empty()
        ch_hsmetrics = channel.empty()
        ch_qualimap  = channel.empty()
        ch_ngsbits   = channel.empty()

        PICARD_COLLECTMULTIPLEMETRICS (ch_bam_bai, ch_genome_fasta, ch_genome_fai)

        ch_bam_bai
            .combine(ch_bait_intervals)
            .combine(ch_target_intervals)
            .set { ch_hsmetrics_in}

        if (val_target_bed) {
            ch_hsmetrics = PICARD_COLLECTHSMETRICS (ch_hsmetrics_in, ch_genome_fasta, ch_genome_fai, [[],[]]).metrics
        }
        if (!skip_qualimap) {
            ch_qualimap = QUALIMAP_BAMQC (ch_bam, []).results
            ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions)
        }

        TIDDIT_COV (ch_bam, [[],[]]) // 2nd pos. arg is req. only for cram input

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
                ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS_WG.out.versions, PICARD_COLLECTWGSMETRICS_Y.out.versions)
            } else if (val_aligner.equals("sentieon")) {
                SENTIEON_WGSMETRICS_WG ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_intervals_wgs.map{ interval -> [[:], interval]} )
                SENTIEON_WGSMETRICS_Y ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_intervals_y.map{ interval -> [[:], interval]} )
                ch_cov      = SENTIEON_WGSMETRICS_WG.out.wgs_metrics
                ch_cov_y    = SENTIEON_WGSMETRICS_Y.out.wgs_metrics
                ch_versions = ch_versions.mix(SENTIEON_WGSMETRICS_WG.out.versions, SENTIEON_WGSMETRICS_Y.out.versions)
            }
        }
        // Check sex
        if (!skip_ngsbits) {
            NGSBITS_SAMPLEGENDER(ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_ngsbits_method)
            ch_ngsbits  = NGSBITS_SAMPLEGENDER.out.tsv
            ch_versions = ch_versions.mix(NGSBITS_SAMPLEGENDER.out.versions)
        }
        // Check contamination
        ch_svd_in = ch_svd_ud.combine(ch_svd_mu).combine(ch_svd_bed).collect()
        VERIFYBAMID_VERIFYBAMID2(ch_bam_bai, ch_svd_in, [], ch_genome_fasta.map {_meta, fasta-> fasta})

        ch_versions = ch_versions.mix(CHROMOGRAPH_COV.out.versions)
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)
        ch_versions = ch_versions.mix(TIDDIT_COV.out.versions)
        ch_versions = ch_versions.mix(UCSC_WIGTOBIGWIG.out.versions)
        ch_versions = ch_versions.mix(VERIFYBAMID_VERIFYBAMID2.out.versions)
        ch_versions = ch_versions.mix(SAMBAMBA_DEPTH.out.versions)

    emit:
        multiple_metrics = PICARD_COLLECTMULTIPLEMETRICS.out.metrics // channel: [ val(meta), path(metrics) ]
        hs_metrics       = ch_hsmetrics                              // channel: [ val(meta), path(metrics) ]
        qualimap_results = ch_qualimap                               // channel: [ val(meta), path(qualimap_dir) ]
        tiddit_wig       = TIDDIT_COV.out.wig                        // channel: [ val(meta), path(wig) ]
        bigwig           = UCSC_WIGTOBIGWIG.out.bw                   // channel: [ val(meta), path(bw) ]
        d4               = MOSDEPTH.out.per_base_d4                  // channel: [ val(meta), path(d4) ]
        global_dist      = MOSDEPTH.out.global_txt                   // channel: [ val(meta), path(txt) ]
        sex_check        = ch_ngsbits                                // channel: [ val(meta), path(tsv) ]
        self_sm          = VERIFYBAMID_VERIFYBAMID2.out.self_sm      // channel: [ val(meta), path(selfSM) ]
        cov              = ch_cov                                    // channel: [ val(meta), path(metrics) ]
        cov_y            = ch_cov_y                                  // channel: [ val(meta), path(metrics) ]
        versions         = ch_versions                               // channel: [ path(versions.yml) ]
}
