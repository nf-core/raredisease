//
// Align reads to the mitochondrial genome
//

include { ALIGN_MT                   } from '../align_MT'
include { ALIGN_MT as ALIGN_MT_SHIFT } from '../align_MT'
include { CONVERT_MT_BAM_TO_FASTQ    } from '../convert_mt_bam_to_fastq'

workflow ALIGN_MITOCHONDRIA {
    take:
        ch_genome_marked_bam_bai // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_dictionary     // channel: [mandatory] [ val(meta), path(dict) ]
        ch_genome_fai            // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_fasta          // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_mt_bwaindex           // channel: [mandatory] [ val(meta), path(index) ]
        ch_mt_bwamem2index       // channel: [mandatory] [ val(meta), path(index) ]
        ch_mt_dictionary         // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mt_fai                // channel: [mandatory] [ val(meta), path(fai) ]
        ch_mt_fasta              // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_mtshift_bwaindex      // channel: [mandatory] [ val(meta), path(index) ]
        ch_mtshift_bwamem2index  // channel: [mandatory] [ val(meta), path(index) ]
        ch_mtshift_dictionary    // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mtshift_fai           // channel: [mandatory] [ val(meta), path(fai) ]
        ch_mtshift_fasta         // channel: [mandatory] [ val(meta), path(fasta) ]
        val_mt_aligner           //  string:  'bwa', 'bwamem2', or 'sentieon'

    main:
        CONVERT_MT_BAM_TO_FASTQ (
            ch_genome_marked_bam_bai,
            ch_genome_dictionary,
            ch_genome_fai,
            ch_genome_fasta
        )
        ch_mt_fastq = CONVERT_MT_BAM_TO_FASTQ.out.fastq

        ALIGN_MT (
            ch_mt_bwaindex,
            ch_mt_bwamem2index,
            ch_mt_dictionary,
            ch_mt_fai,
            ch_mt_fasta,
            ch_mt_fastq,
            CONVERT_MT_BAM_TO_FASTQ.out.ubam,
            val_mt_aligner
        )

        ALIGN_MT_SHIFT (
            ch_mtshift_bwaindex,
            ch_mtshift_bwamem2index,
            ch_mtshift_dictionary,
            ch_mtshift_fai,
            ch_mtshift_fasta,
            ch_mt_fastq,
            CONVERT_MT_BAM_TO_FASTQ.out.ubam,
            val_mt_aligner
        )

        ch_mt_bam_bai                = CONVERT_MT_BAM_TO_FASTQ.out.bam_bai // Used for subsampling and SV calling
        ch_mt_bam_bai_gatksubwf      = ALIGN_MT.out.marked_bam
                                        .join(ALIGN_MT.out.marked_bai, failOnMismatch:true, failOnDuplicate:true) // Only for SNV calling
        ch_mtshift_bam_bai_gatksubwf = ALIGN_MT_SHIFT.out.marked_bam
                                        .join(ALIGN_MT_SHIFT.out.marked_bai, failOnMismatch:true, failOnDuplicate:true) // Only for SNV calling

    emit:
        mt_fastq                  = ch_mt_fastq                  // channel: [ val(meta), [ path(fastq) ] ]
        mt_bam_bai                = ch_mt_bam_bai                // channel: [ val(meta), path(bam), path(bai) ]
        mt_bam_bai_gatksubwf      = ch_mt_bam_bai_gatksubwf      // channel: [ val(meta), path(bam), path(bai) ]
        mtshift_bam_bai_gatksubwf = ch_mtshift_bam_bai_gatksubwf // channel: [ val(meta), path(bam), path(bai) ]
}
