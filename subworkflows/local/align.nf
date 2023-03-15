//
// Map to reference
//

include { ALIGN_BWAMEM2  } from './alignment/align_bwamem2'
include { ALIGN_SENTIEON } from './alignment/align_sentieon'

workflow ALIGN {
    take:
        ch_reads_input     // channel: [mandatory] [ val(meta), [path(reads)]  ]
        ch_fasta           // channel: [mandatory] [ path(fasta) ]
        ch_fai             // channel: [mandatory] [ path(fai) ]
        ch_index_bwa       // channel: [mandatory] [ val(meta), path(index) ]
        ch_index_bwamem2   // channel: [mandatory] [ val(meta), path(index) ]
        ch_known_dbsnp     // channel: [optional; used by sentieon] [ path(known_dbsnp) ]
        ch_known_dbsnp_tbi // channel: [optional; used by sentieon] [ path(known_dbsnp_tbi) ]
        val_platform       // string:  [mandatory] illumina or a different technology

    main:
        ch_versions   = Channel.empty()

        ALIGN_BWAMEM2 (
            ch_reads_input,
            ch_index_bwamem2,
            ch_fasta,
            ch_fai,
            val_platform
        )

        ALIGN_SENTIEON (
            ch_reads_input,
            ch_fasta,
            ch_fai,
            ch_index_bwa,
            ch_known_dbsnp,
            ch_known_dbsnp_tbi,
            val_platform
        )

        ch_marked_bam = Channel.empty().mix(ALIGN_BWAMEM2.out.marked_bam, ALIGN_SENTIEON.out.marked_bam)
        ch_marked_bai = Channel.empty().mix(ALIGN_BWAMEM2.out.marked_bai, ALIGN_SENTIEON.out.marked_bai)
        ch_bam_bai    = ch_marked_bam.join(ch_marked_bai, by: [0])
        ch_versions   = Channel.empty().mix(ALIGN_BWAMEM2.out.versions, ALIGN_SENTIEON.out.versions)

    emit:
        marked_bam = ch_marked_bam // channel: [ val(meta), path(bam) ]
        marked_bai = ch_marked_bai // channel: [ val(meta), path(bai) ]
        bam_bai    = ch_bam_bai    // channel: [ val(meta), path(bam), path(bai) ]
        versions   = ch_versions   // channel: [ path(versions.yml) ]
}
