//
// A nested subworkflow to call structural variants.
//

// CHANGE: swap this example line for the real subworkflow
params.bwamem2_idx_options = [:]

// CHANGE: swap this example line for the real subworkflow
include { PREPARE_GENOME } from './prepare_genome' addParams(
    options: params.bwamem2_idx_options
)

workflow CALL_STRUCTURAL_VARIANTS {

    take:
        bam     // channel: [ val(meta), path(bam) ]
        fasta   // channel: [ path(genome.fasta) ]
        fai   // channel: [ path(genome.fai) ]

    main:
        ch_versions = Channel.empty()

    emit:
        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
