//
// Map to reference
//

include { ALIGN_BWAMEM2  } from '../nf-core/align_bwamem2'
include { ALIGN_SENTIEON } from './align_sentieon'

workflow ALIGN {
    take:
        aligner         // string:  params.aligner
        reads_input     // channel: [ val(meta), reads_input  ]
        fasta           // channel: [genome.fasta]
        fai             // channel: [genome.fai]
        index_bwa       // channel: [ /path/to/bwamem2/index/  ]
        index_bwamem2   // channel: [ /path/to/bwamem2/index/  ]
        known_dbsnp     // channel: [ /path/to/known_dbsnp     ]
        known_dbsnp_tbi // channel: [ /path/to/known_dbsnp_tbi ]

    main:
        ch_versions   = Channel.empty()

        if( aligner == "bwamem2" ) {
            ALIGN_BWAMEM2 ( reads_input, index_bwamem2, fasta, fai )
            ch_marked_bam = ALIGN_BWAMEM2.out.marked_bam
            ch_marked_bai = ALIGN_BWAMEM2.out.marked_bai
            ch_versions = ch_versions.mix(ALIGN_BWAMEM2.out.versions)
        } else if( aligner == "sentieon" ) {
            ALIGN_SENTIEON ( reads_input, fasta, fai, index_bwa, known_dbsnp, known_dbsnp_tbi )
            ch_marked_bam = ALIGN_SENTIEON.out.marked_bam
            ch_marked_bai = ALIGN_SENTIEON.out.marked_bai
            ch_versions = ch_versions.mix(ALIGN_SENTIEON.out.versions)
        } else {
            exit 1, 'Please provide a valid aligner!'
        }

        ch_bam_bai  = ch_marked_bam.join(ch_marked_bai, by: [0])

    emit:
        marked_bam             = ch_marked_bam             // channel: [ val(meta), [ marked_bam ] ]
        marked_bai             = ch_marked_bai             // channel: [ val(meta), [ marked_bai ] ]
        bam_bai                = ch_bam_bai                // channel: [ val(meta), [ marked_bam, marked_bai ] ]

        versions               = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
