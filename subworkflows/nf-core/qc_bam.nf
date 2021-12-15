//
// A quality check subworkflow for processed bams.
//

include { PICARD_COLLECTMULTIPLEMETRICS } from '../../modules/nf-core/modules/picard/collectmultiplemetrics/main'
include { QUALIMAP_BAMQC } from '../../modules/nf-core/modules/qualimap/bamqc/main'
include { TIDDIT_COV } from '../../modules/nf-core/modules/tiddit/cov/main'

workflow QC_BAM {

    take:
        bam     // channel: [ val(meta), path(bam) ]
        fasta   // path: genome.fasta

    main:
        ch_versions = Channel.empty()

        // COLLECT MULTIPLE METRICS
        PICARD_COLLECTMULTIPLEMETRICS ( bam, fasta )
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)

        // QUALIMAP BAMQC
        gff = []
        use_gff = false
        QUALIMAP_BAMQC ( bam, gff, use_gff )
        ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions)

        // TIDDIT COVERAGE
        TIDDIT_COV ( bam, [] ) // 2nd pos. arg is req. only for cram input
        ch_versions = ch_versions.mix(TIDDIT_COV.out.versions)


    emit:
        multiple_metrics        = PICARD_COLLECTMULTIPLEMETRICS.out.metrics     // channel: [ val(meta), path(metrics) ]
        qualimap_results        = QUALIMAP_BAMQC.out.results                    // channel: [ val(meta), path(qualimap files) ]

        versions                = ch_versions.ifEmpty(null)                     // channel: [ versions.yml ]
}
