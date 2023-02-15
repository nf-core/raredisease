//
// A quality check subworkflow for processed bams.
//

include { PICARD_COLLECTMULTIPLEMETRICS } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTHSMETRICS       } from '../../modules/nf-core/picard/collecthsmetrics/main'
include { QUALIMAP_BAMQC                } from '../../modules/nf-core/qualimap/bamqc/main'
include { TIDDIT_COV                    } from '../../modules/nf-core/tiddit/cov/main'
include { MOSDEPTH                      } from '../../modules/nf-core/mosdepth/main'
include { UCSC_WIGTOBIGWIG              } from '../../modules/nf-core/ucsc/wigtobigwig/main'

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

    main:
        ch_versions = Channel.empty()

        PICARD_COLLECTMULTIPLEMETRICS (bam_bai, fasta, fai)

        PICARD_COLLECTHSMETRICS (bam_bai, fasta, fai, bait_intervals, target_intervals)

        QUALIMAP_BAMQC (bam, [])

        TIDDIT_COV (bam, []) // 2nd pos. arg is req. only for cram input

        UCSC_WIGTOBIGWIG (TIDDIT_COV.out.wig, chrom_sizes)

        MOSDEPTH (bam_bai, Channel.value([[], []]), Channel.value([[], []]))

        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
        ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions.first())
        ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())
        ch_versions = ch_versions.mix(TIDDIT_COV.out.versions.first())
        ch_versions = ch_versions.mix(UCSC_WIGTOBIGWIG.out.versions.first())
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

    emit:
        multiple_metrics  = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
        hs_metrics        = PICARD_COLLECTHSMETRICS.out.metrics
        qualimap_results  = QUALIMAP_BAMQC.out.results
        tiddit_wig        = TIDDIT_COV.out.wig
        bigwig            = UCSC_WIGTOBIGWIG.out.bw
        d4                = MOSDEPTH.out.per_base_d4
        versions          = ch_versions.ifEmpty(null)
}
