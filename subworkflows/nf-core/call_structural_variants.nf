//
// A nested subworkflow to call structural variants.
//

workflow CALL_STRUCTURAL_VARIANTS {

    take:
        bam     // channel: [ val(meta), path(bam) ]
        fasta   // channel: [ path(genome.fasta) ]
        fai   // channel: [ path(genome.fai) ]

    main:
        ch_versions = Channel.empty()
}
