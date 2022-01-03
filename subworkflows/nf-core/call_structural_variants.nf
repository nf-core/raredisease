//
// A nested subworkflow to call structural variants.
//

// CHANGE: swap this example line for the real subworkflow
include { CALL_SV_MANTA } from '../local/call_sv_manta'

workflow CALL_STRUCTURAL_VARIANTS {

    take:
        bam     // channel: [ val(meta), path(bam) ]
        fasta   // channel: [ path(genome.fasta) ]
        fai     // channel: [ path(genome.fai) ]

    main:
        ch_versions = Channel.empty()

    emit:
        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
