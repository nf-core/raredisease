//
// Map to reference
//

include { ALIGN_BWAMEM2 } from './align_bwamem2'

workflow ALIGN {
    take:
        aligner     // string: params.aligner
        reads_input // channel: [ val(meta), reads_input ]
        index       // channel: [ /path/to/bwamem2/index/ ]

    main:
        ch_versions = Channel.empty()

        //bwamem2
        ALIGN_BWAMEM2 ( reads_input, index )

        if( aligner == "bwamem2" ) {
            ch_marked_bam = ALIGN_BWAMEM2.out.marked_bam
            ch_marked_bai = ALIGN_BWAMEM2.out.marked_bai
        } else {
            exit 1, 'Please provide a valid aligner!'
        }

        ch_bam_bai  = ch_marked_bam.join(ch_marked_bai, by: [0])
        ch_versions = ch_versions.mix(ALIGN_BWAMEM2.out.versions)

    emit:
        marked_bam             = ch_marked_bam          // channel: [ val(meta), [ marked_bam ] ]
        marked_bai             = ch_marked_bai          // channel: [ val(meta), [ marked_bai ] ]
        bam_bai                = ch_bam_bai             // channel: [ val(meta), [ marked_bam, marked_bai ] ]

        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
