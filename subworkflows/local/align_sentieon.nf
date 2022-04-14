//
// A subworkflow to annotate structural variants.
//

include { SENTIEON_BWAMEM } from '../../modules/local/sentieon/bwamem'

workflow ALIGN_SENTIEON {
    take:
        reads_input // channel: [ val(meta), reads_input ]
        fasta       // path: genome.fasta
        fai         // path: genome.fai
        index       // channel: [ /path/to/bwamem2/index/ ]

    main:
        ch_versions = Channel.empty()

        SENTIEON_BWAMEM ( reads_input, fasta, fai, index )
        ch_versions = ch_versions.mix(SENTIEON_BWAMEM.out.versions)

    emit:
        bam                    = SENTIEON_BWAMEM.out.bam
        bai                    = SENTIEON_BWAMEM.out.bai
        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
