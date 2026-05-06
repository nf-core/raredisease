//
// Align MT
//

include { BWAMEM2_MEM as BWAMEM2_MEM_MT                                     } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWA_MEM as BWA_MEM_MT                                             } from '../../../modules/nf-core/bwa/mem/main'
include { GATK4_MERGEBAMALIGNMENT as GATK4_MERGEBAMALIGNMENT_MT             } from '../../../modules/nf-core/gatk4/mergebamalignment/main'
include { PICARD_ADDORREPLACEREADGROUPS as PICARD_ADDORREPLACEREADGROUPS_MT } from '../../../modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_MARKDUPLICATES as PICARD_MARKDUPLICATES_MT                 } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_MT                                 } from '../../../modules/nf-core/samtools/sort/main'
include { SENTIEON_BWAMEM as SENTIEON_BWAMEM_MT                             } from '../../../modules/nf-core/sentieon/bwamem/main'

workflow ALIGN_MT {
    take:
        ch_bwaindex      // channel: [mandatory for sentieon] [ val(meta), path(index) ]
        ch_bwamem2index  // channel: [mandatory for bwamem2] [ val(meta), path(index) ]
        ch_dict          // channel: [mandatory] [ val(meta), path(dict) ]
        ch_fai           // channel: [mandatory] [ val(meta), path(fai) ]
        ch_fasta         // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_fastq         // channel: [mandatory] [ val(meta), [ path(reads) ] ]
        ch_ubam          // channel: [mandatory] [ val(meta), path(bam) ]
        val_mt_aligner   // string:  'bwa', 'bwamem2', or 'sentieon'

    main:

        if (val_mt_aligner.equals("bwamem2")) {
            BWAMEM2_MEM_MT (ch_fastq, ch_bwamem2index, ch_fasta, true)
            ch_align       = BWAMEM2_MEM_MT.out.bam
        } else if (val_mt_aligner.equals("sentieon")) {
            SENTIEON_BWAMEM_MT ( ch_fastq, ch_bwaindex, ch_fasta, ch_fai )
            ch_align       = SENTIEON_BWAMEM_MT.out.bam_and_bai.map{ meta, bam, _bai -> [meta, bam] }
        } else if (val_mt_aligner.equals("bwa")) {
            BWA_MEM_MT ( ch_fastq, ch_bwaindex, ch_fasta, true )
            ch_align       = BWA_MEM_MT.out.bam
        }
        ch_align
            .join(ch_ubam, failOnMismatch:true, failOnDuplicate:true)
            .set {ch_bam_ubam}

        GATK4_MERGEBAMALIGNMENT_MT (ch_bam_ubam, ch_fasta, ch_dict)

        PICARD_ADDORREPLACEREADGROUPS_MT (GATK4_MERGEBAMALIGNMENT_MT.out.bam, [[:],[]], [[:],[]])

        PICARD_MARKDUPLICATES_MT (PICARD_ADDORREPLACEREADGROUPS_MT.out.bam, ch_fasta, ch_fai)

        SAMTOOLS_SORT_MT (PICARD_MARKDUPLICATES_MT.out.bam, [[:],[]], 'bai')

    emit:
        marked_bai  = SAMTOOLS_SORT_MT.out.bai   // channel: [ val(meta), path(bai) ]
        marked_bam  = SAMTOOLS_SORT_MT.out.bam   // channel: [ val(meta), path(bam) ]
}
