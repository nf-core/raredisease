//
// Align MT
//

include { BWA_MEM as BWA_MEM_MT                                             } from '../../../modules/nf-core/bwa/mem/main'
include { SENTIEON_BWAMEM as SENTIEON_BWAMEM_MT                             } from '../../../modules/nf-core/sentieon/bwamem/main'
include { BWAMEM2_MEM as BWAMEM2_MEM_MT                                     } from '../../../modules/nf-core/bwamem2/mem/main'
include { GATK4_MERGEBAMALIGNMENT as GATK4_MERGEBAMALIGNMENT_MT             } from '../../../modules/nf-core/gatk4/mergebamalignment/main'
include { PICARD_ADDORREPLACEREADGROUPS as PICARD_ADDORREPLACEREADGROUPS_MT } from '../../../modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_MARKDUPLICATES as PICARD_MARKDUPLICATES_MT                 } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MT                               } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_MT                                 } from '../../../modules/nf-core/samtools/sort/main'

workflow ALIGN_MT {
    take:
        ch_fastq        // channel: [mandatory] [ val(meta), [ path(reads) ] ]
        ch_ubam         // channel: [mandatory] [ val(meta), path(bam) ]
        ch_bwaindex     // channel: [mandatory for sentieon] [ val(meta), path(index) ]
        ch_bwamem2index // channel: [mandatory for bwamem2] [ val(meta), path(index) ]
        ch_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_dict         // channel: [mandatory] [ val(meta), path(dict) ]
        ch_fai          // channel: [mandatory] [ val(meta), path(fai) ]

    main:
        ch_versions     = Channel.empty()
        ch_bwa_bam      = Channel.empty()
        ch_bwamem2_bam  = Channel.empty()
        ch_sentieon_bam = Channel.empty()

        if (params.aligner.equals("bwamem2")) {
            BWAMEM2_MEM_MT (ch_fastq, ch_bwamem2index, true)
            ch_bwamem2_bam = BWAMEM2_MEM_MT.out.bam
            ch_versions    = ch_versions.mix(BWAMEM2_MEM_MT.out.versions.first())
        } else if (params.aligner.equals("sentieon")) {
            SENTIEON_BWAMEM_MT ( ch_fastq, ch_bwaindex, ch_fasta, ch_fai )
            ch_sentieon_bam = SENTIEON_BWAMEM_MT.out.bam_and_bai.map{ meta, bam, bai -> [meta, bam] }
            ch_versions     = ch_versions.mix(SENTIEON_BWAMEM_MT.out.versions.first())
        } else if (params.aligner.equals("bwa")) {
            BWA_MEM_MT ( ch_fastq, ch_bwaindex, true )
            ch_bwa_bam      = BWA_MEM_MT.out.bam
            ch_versions     = ch_versions.mix(BWA_MEM_MT.out.versions.first())
        }
        Channel.empty()
            .mix(ch_bwamem2_bam, ch_sentieon_bam, ch_bwa_bam)
            .join(ch_ubam, failOnMismatch:true, failOnDuplicate:true)
            .set {ch_bam_ubam}

        GATK4_MERGEBAMALIGNMENT_MT (ch_bam_ubam, ch_fasta, ch_dict)

        PICARD_ADDORREPLACEREADGROUPS_MT (GATK4_MERGEBAMALIGNMENT_MT.out.bam)

        PICARD_MARKDUPLICATES_MT (PICARD_ADDORREPLACEREADGROUPS_MT.out.bam, ch_fasta, ch_fai)

        SAMTOOLS_SORT_MT (PICARD_MARKDUPLICATES_MT.out.bam)

        SAMTOOLS_INDEX_MT(SAMTOOLS_SORT_MT.out.bam)

        ch_versions = ch_versions.mix(GATK4_MERGEBAMALIGNMENT_MT.out.versions.first())
        ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS_MT.out.versions.first())
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES_MT.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_SORT_MT.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MT.out.versions.first())

    emit:
        marked_bam  = SAMTOOLS_SORT_MT.out.bam   // channel: [ val(meta), path(bam) ]
        marked_bai  = SAMTOOLS_INDEX_MT.out.bai  // channel: [ val(meta), path(bai) ]
        versions    = ch_versions                // channel: [ path(versions.yml) ]
}
