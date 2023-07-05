//
// Map to reference
//

include { ALIGN_BWAMEM2  } from './alignment/align_bwamem2'
include { ALIGN_SENTIEON } from './alignment/align_sentieon'
include { SAMTOOLS_VIEW  } from '../../modules/nf-core/samtools/view/main'

workflow ALIGN {
    take:
        ch_reads_input     // channel: [mandatory] [ val(meta), [path(reads)]  ]
        ch_genome_fasta    // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai      // channel: [mandatory] [ val(meta), path(fai) ]
        ch_bwa_index       // channel: [mandatory] [ val(meta), path(index) ]
        ch_bwamem2_index   // channel: [mandatory] [ val(meta), path(index) ]
        ch_known_dbsnp     // channel: [optional; used by sentieon] [ path(known_dbsnp) ]
        ch_known_dbsnp_tbi // channel: [optional; used by sentieon] [ path(known_dbsnp_tbi) ]
        val_platform       // string:  [mandatory] illumina or a different technology

    main:
        ch_versions   = Channel.empty()

        ALIGN_BWAMEM2 (             // Triggered when params.aligner is set as bwamem2
            ch_reads_input,
            ch_bwamem2_index,
            ch_genome_fasta,
            ch_genome_fai,
            val_platform
        )

        ALIGN_SENTIEON (            // Triggered when params.aligner is set as sentieon
            ch_reads_input,
            ch_genome_fasta,
            ch_genome_fai,
            ch_bwa_index,
            ch_known_dbsnp,
            ch_known_dbsnp_tbi,
            val_platform
        )

        ch_marked_bam = Channel.empty().mix(ALIGN_BWAMEM2.out.marked_bam, ALIGN_SENTIEON.out.marked_bam)
        ch_marked_bai = Channel.empty().mix(ALIGN_BWAMEM2.out.marked_bai, ALIGN_SENTIEON.out.marked_bai)
        ch_bam_bai    = ch_marked_bam.join(ch_marked_bai, failOnMismatch:true, failOnDuplicate:true)

        SAMTOOLS_VIEW( ch_bam_bai, ch_genome_fasta, [] )

        ch_versions   = Channel.empty().mix(ALIGN_BWAMEM2.out.versions, ALIGN_SENTIEON.out.versions)

    emit:
        marked_bam = ch_marked_bam // channel: [ val(meta), path(bam) ]
        marked_bai = ch_marked_bai // channel: [ val(meta), path(bai) ]
        bam_bai    = ch_bam_bai    // channel: [ val(meta), path(bam), path(bai) ]
        versions   = ch_versions   // channel: [ path(versions.yml) ]
}
