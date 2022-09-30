//
// A quality check subworkflow for processed bams.
//

include { PICARD_COLLECTMULTIPLEMETRICS } from '../../modules/nf-core/modules/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTHSMETRICS       } from '../../modules/nf-core/modules/picard/collecthsmetrics/main'
include { QUALIMAP_BAMQC                } from '../../modules/nf-core/modules/qualimap/bamqc/main'
include { TIDDIT_COV                    } from '../../modules/nf-core/modules/tiddit/cov/main'
include { MOSDEPTH                      } from '../../modules/nf-core/modules/mosdepth/main'
include { UCSC_WIGTOBIGWIG              } from '../../modules/nf-core/modules/ucsc/wigtobigwig/main'

workflow QC_BAM {

    take:
        bam              // channel: [ val(meta), path(bam) ]
        bai              // channel: [ val(meta), path(bai) ]
        fasta            // path: genome.fasta
        fai              // path: genome.fasta.fai
        bait_intervals   // path: bait.intervals_list
        target_intervals // path: target.intervals_list
        chrom_sizes      // path: chrom.sizes

    main:
        ch_versions = Channel.empty()

        // COLLECT MULTIPLE METRICS
        PICARD_COLLECTMULTIPLEMETRICS ( bam, fasta, fai )
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)

        // COLLECT HS METRICS
        PICARD_COLLECTHSMETRICS ( bam, fasta, fai, bait_intervals, target_intervals )
        ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions)

        // QUALIMAP BAMQC
        gff = []
        QUALIMAP_BAMQC ( bam, gff )
        ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions)

        // TIDDIT COVERAGE
        TIDDIT_COV ( bam, [] ) // 2nd pos. arg is req. only for cram input
        UCSC_WIGTOBIGWIG ( TIDDIT_COV.out.wig, chrom_sizes )
        ch_versions = ch_versions.mix(TIDDIT_COV.out.versions)

        // MOSDEPTH
        mosdepth_input_bams = bam.join(bai, by: [0])
        MOSDEPTH (mosdepth_input_bams,[],[])
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)

    emit:
        multiple_metrics        = PICARD_COLLECTMULTIPLEMETRICS.out.metrics     // channel: [ val(meta), path(metrics) ]
        hs_metrics              = PICARD_COLLECTHSMETRICS.out.metrics           // channel: [ val(meta), path(metrics) ]
        qualimap_results        = QUALIMAP_BAMQC.out.results                    // channel: [ val(meta), path(qualimap files) ]
        tiddit_wig              = TIDDIT_COV.out.wig                            // channel: [ val(meta), path(*.wig) ]
        bigwig                  = UCSC_WIGTOBIGWIG.out.bw                       // channel: [ val(meta), path(*.bw) ]
        d4                      = MOSDEPTH.out.per_base_d4                      // channel: [ val(meta), path(*.d4) ]

        versions                = ch_versions.ifEmpty(null)                     // channel: [ versions.yml ]
}
