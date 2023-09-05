//
// Map to reference
//

include { ALIGN_BWAMEM2              } from './alignment/align_bwamem2'
include { ALIGN_SENTIEON             } from './alignment/align_sentieon'
include { SAMTOOLS_VIEW              } from '../../modules/nf-core/samtools/view/main'
include { ALIGN_MT                   } from './alignment/align_MT'
include { ALIGN_MT as ALIGN_MT_SHIFT } from './alignment/align_MT'
include { CONVERT_MT_BAM_TO_FASTQ    } from './convert_mt_bam_to_fastq'

workflow ALIGN {
    take:
        ch_reads                 // channel: [mandatory] [ val(meta), [path(reads)]  ]
        ch_genome_fasta          // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai            // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_bwaindex       // channel: [mandatory] [ val(meta), path(index) ]
        ch_genome_bwamem2index   // channel: [mandatory] [ val(meta), path(index) ]
        ch_genome_dictionary     // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mtshift_bwaindex      // channel: [mandatory] [ val(meta), path(index) ]
        ch_mtshift_bwamem2index  // channel: [mandatory] [ val(meta), path(index) ]
        ch_mtshift_fasta         // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_mtshift_dictionary    // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mtshift_fai           // channel: [mandatory] [ val(meta), path(fai) ]
        val_platform             // string:  [mandatory] illumina or a different technology

    main:
        ch_versions   = Channel.empty()

        ALIGN_BWAMEM2 (             // Triggered when params.aligner is set as bwamem2
            ch_reads,
            ch_genome_bwamem2index,
            ch_genome_fasta,
            ch_genome_fai,
            val_platform
        )

        ALIGN_SENTIEON (            // Triggered when params.aligner is set as sentieon
            ch_reads,
            ch_genome_fasta,
            ch_genome_fai,
            ch_genome_bwaindex,
            val_platform
        )

        ch_genome_marked_bam = Channel.empty().mix(ALIGN_BWAMEM2.out.marked_bam, ALIGN_SENTIEON.out.marked_bam)
        ch_genome_marked_bai = Channel.empty().mix(ALIGN_BWAMEM2.out.marked_bai, ALIGN_SENTIEON.out.marked_bai)
        ch_genome_bam_bai    = ch_genome_marked_bam.join(ch_genome_marked_bai, failOnMismatch:true, failOnDuplicate:true)

        // PREPARING READS FOR MT ALIGNMENT
        CONVERT_MT_BAM_TO_FASTQ (
            ch_genome_bam_bai,
            ch_genome_fasta,
            ch_genome_fai,
            ch_genome_dictionary
        )

        ALIGN_MT (
            CONVERT_MT_BAM_TO_FASTQ.out.fastq,
            CONVERT_MT_BAM_TO_FASTQ.out.bam,
            ch_genome_bwaindex,
            ch_genome_bwamem2index,
            ch_genome_fasta,
            ch_genome_dictionary,
            ch_genome_fai
        )

        ALIGN_MT_SHIFT (
            CONVERT_MT_BAM_TO_FASTQ.out.fastq,
            CONVERT_MT_BAM_TO_FASTQ.out.bam,
            ch_mtshift_bwaindex,
            ch_mtshift_bwamem2index,
            ch_mtshift_fasta,
            ch_mtshift_dictionary,
            ch_mtshift_fai
        )

        ch_mt_marked_bam = ALIGN_MT.out.marked_bam
        ch_mt_marked_bai = ALIGN_MT.out.marked_bai
        ch_mt_bam_bai    = ch_mt_marked_bam.join(ch_mt_marked_bai, failOnMismatch:true, failOnDuplicate:true)

        ch_mtshift_marked_bam = ALIGN_MT_SHIFT.out.marked_bam
        ch_mtshift_marked_bai = ALIGN_MT_SHIFT.out.marked_bai
        ch_mtshift_bam_bai    = ch_mtshift_marked_bam.join(ch_mtshift_marked_bai, failOnMismatch:true, failOnDuplicate:true)

        SAMTOOLS_VIEW( ch_genome_bam_bai, ch_genome_fasta, [] )

        ch_versions   = Channel.empty().mix(ALIGN_BWAMEM2.out.versions,
                                            ALIGN_SENTIEON.out.versions,
                                            ALIGN_MT.out.versions,
                                            ALIGN_MT_SHIFT.out.versions,
                                            CONVERT_MT_BAM_TO_FASTQ.out.versions,
                                            SAMTOOLS_VIEW.out.versions)

    emit:
        genome_marked_bam  = ch_genome_marked_bam  // channel: [ val(meta), path(bam) ]
        genome_marked_bai  = ch_genome_marked_bai  // channel: [ val(meta), path(bai) ]
        genome_bam_bai     = ch_genome_bam_bai     // channel: [ val(meta), path(bam), path(bai) ]
        mt_marked_bam      = ch_mt_marked_bam      // channel: [ val(meta), path(bam) ]
        mt_marked_bai      = ch_mt_marked_bai      // channel: [ val(meta), path(bai) ]
        mt_bam_bai         = ch_mt_bam_bai         // channel: [ val(meta), path(bam), path(bai) ]
        mtshift_marked_bam = ch_mtshift_marked_bam // channel: [ val(meta), path(bam) ]
        mtshift_marked_bai = ch_mtshift_marked_bai // channel: [ val(meta), path(bai) ]
        mtshift_bam_bai    = ch_mtshift_bam_bai    // channel: [ val(meta), path(bam), path(bai) ]
        versions           = ch_versions           // channel: [ path(versions.yml) ]
}
