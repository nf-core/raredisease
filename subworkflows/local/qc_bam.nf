//
// A quality check subworkflow for processed bams.
//

include { PICARD_COLLECTMULTIPLEMETRICS                          } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTHSMETRICS                                } from '../../modules/nf-core/picard/collecthsmetrics/main'
include { QUALIMAP_BAMQC                                         } from '../../modules/nf-core/qualimap/bamqc/main'
include { TIDDIT_COV                                             } from '../../modules/nf-core/tiddit/cov/main'
include { MOSDEPTH                                               } from '../../modules/nf-core/mosdepth/main'
include { UCSC_WIGTOBIGWIG                                       } from '../../modules/nf-core/ucsc/wigtobigwig/main'
include { PICARD_COLLECTWGSMETRICS as PICARD_COLLECTWGSMETRICS   } from '../../modules/nf-core/picard/collectwgsmetrics/main'
include { PICARD_COLLECTWGSMETRICS as PICARD_COLLECTWGSMETRICS_Y } from '../../modules/nf-core/picard/collectwgsmetrics/main'
include { SENTIEON_WGSMETRICSALGO as SENTIEON_WGSMETRICS         } from '../../modules/local/sentieon/wgsmetricsalgo'
include { SENTIEON_WGSMETRICSALGO as SENTIEON_WGSMETRICS_Y       } from '../../modules/local/sentieon/wgsmetricsalgo'

workflow QC_BAM {

    take:
        bam              // channel: [ val(meta), path(bam) ]
        bai              // channel: [ val(meta), path(bai) ]
        bam_bai
        fasta            // path: genome.fasta
        fai              // path: genome.fasta.fai
        bait_intervals   // path: bait.intervals_list
        target_intervals // path: target.intervals_list
        chrom_sizes      // path: chrom.sizes
        intervals_wgs    // path: genome.intervals_wgs
        intervals_y      // path: genome.intervals_y
        aligner          // string: params.aligner

    main:
        ch_versions = Channel.empty()

        PICARD_COLLECTMULTIPLEMETRICS (bam_bai, fasta, fai)

        PICARD_COLLECTHSMETRICS (bam_bai, fasta, fai, bait_intervals, target_intervals)

        QUALIMAP_BAMQC (bam, [])

        TIDDIT_COV (bam, []) // 2nd pos. arg is req. only for cram input

        UCSC_WIGTOBIGWIG (TIDDIT_COV.out.wig, chrom_sizes)

        MOSDEPTH (bam_bai, [], [])

        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
        ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions.first())
        ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())
        ch_versions = ch_versions.mix(TIDDIT_COV.out.versions.first())
        ch_versions = ch_versions.mix(UCSC_WIGTOBIGWIG.out.versions.first())
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

        // COLLECT WGS METRICS
        PICARD_COLLECTWGSMETRICS ( bam_bai, fasta, fai, intervals_wgs )
        PICARD_COLLECTWGSMETRICS_Y ( bam_bai, fasta, fai, intervals_y )
        
        SENTIEON_WGSMETRICS ( bam_bai, fasta, fai)
        SENTIEON_WGSMETRICS_Y ( bam_bai, fasta, fai)

        ch_cov   = Channel.empty().mix(PICARD_COLLECTWGSMETRICS.out.metrics, SENTIEON_WGSMETRICS.out.wgs_metrics)
        ch_cov_y = Channel.empty().mix(PICARD_COLLECTWGSMETRICS_Y.out.metrics, SENTIEON_WGSMETRICS_Y.out.wgs_metrics)
        
        ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS.out.versions, SENTIEON_WGSMETRICS.out.versions)
        ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS_Y.out.versions, SENTIEON_WGSMETRICS_Y.out.versions)

    emit:
        multiple_metrics        = PICARD_COLLECTMULTIPLEMETRICS.out.metrics     // channel: [ val(meta), path(metrics) ]
        hs_metrics              = PICARD_COLLECTHSMETRICS.out.metrics           // channel: [ val(meta), path(metrics) ]
        qualimap_results        = QUALIMAP_BAMQC.out.results                    // channel: [ val(meta), path(qualimap files) ]
        tiddit_wig              = TIDDIT_COV.out.wig                            // channel: [ val(meta), path(*.wig) ]
        bigwig                  = UCSC_WIGTOBIGWIG.out.bw                       // channel: [ val(meta), path(*.bw) ]
        d4                      = MOSDEPTH.out.per_base_d4                      // channel: [ val(meta), path(*.d4) ]
        cov                     = ch_cov                                        // channel: [ val(meta), path(metrics) ]
        cov_y                   = ch_cov_y                                      // channel: [ val(meta), path(metrics) ]

        versions                = ch_versions.ifEmpty(null)                     // channel: [ versions.yml ]
}
