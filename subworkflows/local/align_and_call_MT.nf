//
// Allign and call MT
//

include { BWAMEM2_MEM as BWAMEM2_MEM_MT                                     } from '../../modules/nf-core/bwamem2/mem/main'
include { GATK4_MERGEBAMALIGNMENT as GATK4_MERGEBAMALIGNMENT_MT             } from '../../modules/nf-core/gatk4/mergebamalignment/main'
include { PICARD_ADDORREPLACEREADGROUPS as PICARD_ADDORREPLACEREADGROUPS_MT } from '../../modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_MARKDUPLICATES as PICARD_MARKDUPLICATES_MT                 } from '../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MT                               } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_MT                                 } from '../../modules/nf-core/samtools/sort/main'
include { HAPLOCHECK as HAPLOCHECK_MT                                       } from '../../modules/nf-core/haplocheck/main'
include { GATK4_MUTECT2 as GATK4_MUTECT2_MT                                 } from '../../modules/nf-core/gatk4/mutect2/main'



workflow ALIGN_AND_CALL_MT {
    take:
        fastq         // channel: [ val(meta), path('*.fastq.gz') ]
        ubam          // channel: [ val(meta), path('*.bam') ]
        index         // channel: [ /path/to/bwamem2/index/ ]
        fasta         // channel: [ genome.fasta ]
        dict          // channel: [ genome.dict ]
        fai           // channel: [ genome.fai ]
        intervals_mt  // channel: [ file(non_control_region.chrM.interval_list) ]

    main:
        ch_versions = Channel.empty()

        // Outputs bam files
        BWAMEM2_MEM_MT ( fastq , index, true)
        ch_versions    = ch_versions.mix(BWAMEM2_MEM_MT.out.versions.first())
        ch_mt_bam      =  BWAMEM2_MEM_MT.out.bam
        ch_fastq_ubam  = ch_mt_bam.join(ubam, by: [0])

        // Merges bam files
        GATK4_MERGEBAMALIGNMENT_MT (ch_fastq_ubam, fasta, dict )
        ch_versions = ch_versions.mix(GATK4_MERGEBAMALIGNMENT_MT.out.versions.first())

        // Add read group to merged bam file
        PICARD_ADDORREPLACEREADGROUPS_MT ( GATK4_MERGEBAMALIGNMENT_MT.out.bam )
        ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS_MT.out.versions.first())

        // Marks duplicates
        PICARD_MARKDUPLICATES_MT (PICARD_ADDORREPLACEREADGROUPS_MT.out.bam )
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES_MT.out.versions.first())

        // Sort bam file
        SAMTOOLS_SORT_MT (PICARD_MARKDUPLICATES_MT.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_SORT_MT.out.versions.first())

        // Index bam file
        SAMTOOLS_INDEX_MT(SAMTOOLS_SORT_MT.out.bam)
        ch_sort_index_bam=SAMTOOLS_SORT_MT.out.bam.join(SAMTOOLS_INDEX_MT.out.bai, by: [0])
        ch_sort_index_bam_intervals_mt=ch_sort_index_bam.combine(intervals_mt)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MT.out.versions.first())

        // Calls variants with Mutect2
        GATK4_MUTECT2_MT (ch_sort_index_bam_intervals_mt, fasta, fai, dict, [], [], [],[])
        ch_versions = ch_versions.mix(GATK4_MUTECT2_MT.out.versions.first())

        // Haplocheck
        // TODO: probably it will be outside this subworkflow as we want to run
        // with the VCF with the variants from the shifted alignment (to solve the mt circularity issue)
        HAPLOCHECK_MT ( GATK4_MUTECT2_MT.out.vcf )
        ch_versions = ch_versions.mix(HAPLOCHECK_MT.out.versions.first())


    emit:
        vcf      = GATK4_MUTECT2_MT.out.vcf
        tbi      = GATK4_MUTECT2_MT.out.tbi
        txt      = HAPLOCHECK_MT.out.txt
        html     = HAPLOCHECK_MT.out.html
        versions = ch_versions // channel: [ versions.yml ]
}
