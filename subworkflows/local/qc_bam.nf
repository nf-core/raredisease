//
// A quality check subworkflow for processed bams.
//

include { PICARD_COLLECTMULTIPLEMETRICS                          } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTHSMETRICS                                } from '../../modules/nf-core/picard/collecthsmetrics/main'
include { CHROMOGRAPH as CHROMOGRAPH_COV                         } from '../../modules/nf-core/chromograph/main'
include { QUALIMAP_BAMQC                                         } from '../../modules/nf-core/qualimap/bamqc/main'
include { TIDDIT_COV                                             } from '../../modules/nf-core/tiddit/cov/main'
include { MOSDEPTH                                               } from '../../modules/nf-core/mosdepth/main'
include { UCSC_WIGTOBIGWIG                                       } from '../../modules/nf-core/ucsc/wigtobigwig/main'
include { PICARD_COLLECTWGSMETRICS as PICARD_COLLECTWGSMETRICS   } from '../../modules/nf-core/picard/collectwgsmetrics/main'
include { PICARD_COLLECTWGSMETRICS as PICARD_COLLECTWGSMETRICS_Y } from '../../modules/nf-core/picard/collectwgsmetrics/main'
include { SENTIEON_WGSMETRICS                                    } from '../../modules/nf-core/sentieon/wgsmetrics/main'
include { SENTIEON_WGSMETRICS as SENTIEON_WGSMETRICS_Y           } from '../../modules/nf-core/sentieon/wgsmetrics/main'
include { NGSBITS_SAMPLEGENDER                                   } from '../../modules/nf-core/ngsbits/samplegender/main'

workflow QC_BAM {

    take:
        ch_bam                      // channel: [mandatory] [ val(meta), path(bam) ]
        ch_bai                      // channel: [mandatory] [ val(meta), path(bai) ]
        ch_bam_bai                  // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_fasta             // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai               // channel: [mandatory] [ val(meta), path(fai) ]
        ch_bait_intervals           // channel: [mandatory] [ path(intervals_list) ]
        ch_target_intervals         // channel: [mandatory] [ path(intervals_list) ]
        ch_chrom_sizes              // channel: [mandatory] [ path(sizes) ]
        ch_intervals_wgs            // channel: [mandatory] [ path(intervals) ]
        ch_intervals_y              // channel: [mandatory] [ path(intervals) ]
        ngsbits_samplegender_method // channel [val(method)]

    main:
        ch_versions = Channel.empty()
        ch_qualimap = Channel.empty()

        PICARD_COLLECTMULTIPLEMETRICS (ch_bam_bai, ch_genome_fasta, ch_genome_fai)

        ch_bam_bai
            .combine(ch_bait_intervals)
            .combine(ch_target_intervals)
            .set { ch_hsmetrics_in}

        PICARD_COLLECTHSMETRICS (ch_hsmetrics_in, ch_genome_fasta, ch_genome_fai, [[],[]])

        if (!params.skip_qualimap) {
            ch_qualimap = QUALIMAP_BAMQC (ch_bam, []).results
            ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())
        }

        TIDDIT_COV (ch_bam, [[],[]]) // 2nd pos. arg is req. only for cram input

        UCSC_WIGTOBIGWIG (TIDDIT_COV.out.wig, ch_chrom_sizes)

        CHROMOGRAPH_COV([[:],[]], TIDDIT_COV.out.wig, [[:],[]], [[:],[]], [[:],[]], [[:],[]], [[:],[]])

        ch_bam_bai.map{ meta, bam, bai -> [meta, bam, bai, []]}.set{ch_mosdepth_in}
        MOSDEPTH (ch_mosdepth_in, ch_genome_fasta)

        // COLLECT WGS METRICS
        PICARD_COLLECTWGSMETRICS ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_intervals_wgs )
        PICARD_COLLECTWGSMETRICS_Y ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_intervals_y )

        SENTIEON_WGSMETRICS ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_intervals_wgs.map{ interval -> [[:], interval]} )
        SENTIEON_WGSMETRICS_Y ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_intervals_y.map{ interval -> [[:], interval]} )

        // Check sex
        NGSBITS_SAMPLEGENDER(ch_bam_bai, ch_genome_fasta, ch_genome_fai, ngsbits_samplegender_method)

        ch_cov   = Channel.empty().mix(PICARD_COLLECTWGSMETRICS.out.metrics, SENTIEON_WGSMETRICS.out.wgs_metrics)
        ch_cov_y = Channel.empty().mix(PICARD_COLLECTWGSMETRICS_Y.out.metrics, SENTIEON_WGSMETRICS_Y.out.wgs_metrics)

        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
        ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions.first())
        ch_versions = ch_versions.mix(TIDDIT_COV.out.versions.first())
        ch_versions = ch_versions.mix(UCSC_WIGTOBIGWIG.out.versions.first())
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())
        ch_versions = ch_versions.mix(NGSBITS_SAMPLEGENDER.out.versions.first())
        ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS.out.versions.first(), SENTIEON_WGSMETRICS.out.versions.first())
        ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS_Y.out.versions.first(), SENTIEON_WGSMETRICS_Y.out.versions.first())

    emit:
        multiple_metrics = PICARD_COLLECTMULTIPLEMETRICS.out.metrics // channel: [ val(meta), path(metrics) ]
        hs_metrics       = PICARD_COLLECTHSMETRICS.out.metrics       // channel: [ val(meta), path(metrics) ]
        qualimap_results = ch_qualimap                               // channel: [ val(meta), path(qualimap_dir) ]
        tiddit_wig       = TIDDIT_COV.out.wig                        // channel: [ val(meta), path(wig) ]
        bigwig           = UCSC_WIGTOBIGWIG.out.bw                   // channel: [ val(meta), path(bw) ]
        d4               = MOSDEPTH.out.per_base_d4                  // channel: [ val(meta), path(d4) ]
        global_dist      = MOSDEPTH.out.global_txt                   // channel: [ val(meta), path(txt) ]
        sex_check        = NGSBITS_SAMPLEGENDER.out.tsv              // channel: [val(meta), path(tsv) ]
        cov              = ch_cov                                    // channel: [ val(meta), path(metrics) ]
        cov_y            = ch_cov_y                                  // channel: [ val(meta), path(metrics) ]
        versions         = ch_versions                               // channel: [ path(versions.yml) ]
}
