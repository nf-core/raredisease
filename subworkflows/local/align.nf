//
// Map to reference
//

include { ALIGN_BWAMEM2  } from './alignment/align_bwamem2'
include { ALIGN_SENTIEON } from './alignment/align_sentieon'

workflow ALIGN {
    take:
        reads_input     // channel: [ val(meta), reads_input  ]
        fasta           // channel: [genome.fasta]
        fai             // channel: [genome.fai]
        index_bwa       // channel: [ /path/to/bwamem2/index/  ]
        index_bwamem2   // channel: [ /path/to/bwamem2/index/  ]
        known_dbsnp     // channel: [ /path/to/known_dbsnp     ]
        known_dbsnp_tbi // channel: [ /path/to/known_dbsnp_tbi ]
        platform        // string:  params.platform

    main:
        ch_versions   = Channel.empty()

        ALIGN_BWAMEM2 (
            reads_input,
            index_bwamem2,
            fasta,
            fai,
            platform
        )

        ALIGN_SENTIEON (
            reads_input,
            fasta,
            fai,
            index_bwa,
            known_dbsnp,
            known_dbsnp_tbi,
            platform
        )

        ch_marked_bam = Channel.empty().mix(ALIGN_BWAMEM2.out.marked_bam, ALIGN_SENTIEON.out.marked_bam)
        ch_marked_bai = Channel.empty().mix(ALIGN_BWAMEM2.out.marked_bai, ALIGN_SENTIEON.out.marked_bai)
        ch_bam_bai    = ch_marked_bam.join(ch_marked_bai, by: [0])
        ch_versions   = Channel.empty().mix(ALIGN_BWAMEM2.out.versions, ALIGN_SENTIEON.out.versions)

    emit:
        marked_bam    = ch_marked_bam
        marked_bai    = ch_marked_bai
        bam_bai       = ch_bam_bai
        versions      = ch_versions
}
