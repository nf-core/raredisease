//
// A quality check subworkflow for processed bams.
//

// params.picard_collectmultiplemetrics_options = [:]
// params.qualimap_bamqc_options = [:]

include { PICARD_COLLECTMULTIPLEMETRICS } from '../../modules/nf-core/modules/picard/collectmultiplemetrics/main'  //addParams( options: params.picard_collectmultiplemetrics_options )
include { QUALIMAP_BAMQC } from '../../modules/nf-core/modules/qualimap/bamqc/main'  //addParams( options: params.qualimap_bamqc_options )

workflow QC_BAM {

    take:
        bam     // channel: [ val(meta), path(bam) ]
        fasta   // path: genome.fasta

    main:
        ch_versions = Channel.empty()

        PICARD_COLLECTMULTIPLEMETRICS ( bam, fasta )
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)

        gff = []
        use_gff = false
        QUALIMAP_BAMQC ( bam, gff, use_gff )
        ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions)

    emit:
        multiple_metrics        = PICARD_COLLECTMULTIPLEMETRICS.out.metrics     // channel: [ val(meta), path(metrics) ]
        qualimap_results        = QUALIMAP_BAMQC.out.results                    // channel: [ val(meta), path(qualimap files) ]

        versions                = ch_versions.ifEmpty(null)                     // channel: [ versions.yml ]
}
