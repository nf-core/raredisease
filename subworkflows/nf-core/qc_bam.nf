//
// A quality check subworkflow for processed bams.
//

params.picard_collectmultiplemetrics_options = [:]

include { PICARD_COLLECTMULTIPLEMETRICS } from '../../modules/nf-core/modules/picard/collectmultiplemetrics/main'  addParams( options: params.picard_collectmultiplemetrics_options )

workflow QC_BAM {

    take:
        bam   // channel: [ val(meta), path(bam) ]
        fasta // path: genome.fasta

    main:
        ch_versions = Channel.empty()

        PICARD_COLLECTMULTIPLEMETRICS ( bam, fasta )
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)

    emit:
        multiple_metrics        = PICARD_COLLECTMULTIPLEMETRICS.out.metrics     // channel: [ val(meta), path(metrics) ]

        versions                = ch_versions.ifEmpty(null)                     // channel: [ versions.yml ]
}
