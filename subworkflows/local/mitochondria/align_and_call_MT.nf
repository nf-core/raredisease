//
// Align and call MT
//

include { BWAMEM2_MEM as BWAMEM2_MEM_MT                                     } from '../../modules/nf-core/bwamem2/mem/main'
include { GATK4_MERGEBAMALIGNMENT as GATK4_MERGEBAMALIGNMENT_MT             } from '../../modules/nf-core/gatk4/mergebamalignment/main'
include { PICARD_ADDORREPLACEREADGROUPS as PICARD_ADDORREPLACEREADGROUPS_MT } from '../../modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_MARKDUPLICATES as PICARD_MARKDUPLICATES_MT                 } from '../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MT                               } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_MT                                 } from '../../modules/nf-core/samtools/sort/main'
include { HAPLOCHECK as HAPLOCHECK_MT                                       } from '../../modules/nf-core/haplocheck/main'
include { GATK4_MUTECT2 as GATK4_MUTECT2_MT                                 } from '../../modules/nf-core/gatk4/mutect2/main'
include { GATK4_FILTERMUTECTCALLS as  GATK4_FILTERMUTECTCALLS_MT            } from '../../modules/nf-core/gatk4/filtermutectcalls/main'
include { PICARD_RENAMESAMPLEINVCF as PICARD_RENAMESAMPLEINVCF_MT           } from '../../modules/nf-core/picard/renamesampleinvcf/main'
include { TABIX_TABIX as TABIX_TABIX_MT                                     } from '../../modules/nf-core/tabix/tabix/main'

workflow ALIGN_AND_CALL_MT {
    take:
        fastq         // channel: [ val(meta), path('*.fastq.gz') ]
        ubam          // channel: [ val(meta), path('*.bam') ]
        genome_index  // channel: [ /path/to/bwamem2/index/ ]
        genome_fasta  // channel: [ genome.fasta ]
        genome_dict   // channel: [ genome.dict ]
        genome_fai    // channel: [ genome.fai ]
        intervals_mt  // channel: [ file(non_control_region.chrM.interval_list) ]

    main:
        ch_versions = Channel.empty()

        // Outputs bam files
        BWAMEM2_MEM_MT ( fastq , genome_index, true)
        ch_versions    = ch_versions.mix(BWAMEM2_MEM_MT.out.versions.first())
        ch_mt_bam      =  BWAMEM2_MEM_MT.out.bam
        ch_fastq_ubam  = ch_mt_bam.join(ubam, by: [0])

        // Merges bam files
        GATK4_MERGEBAMALIGNMENT_MT (ch_fastq_ubam, genome_fasta, genome_dict )
        ch_versions = ch_versions.mix(GATK4_MERGEBAMALIGNMENT_MT.out.versions.first())

        // Add read group to merged bam file
        PICARD_ADDORREPLACEREADGROUPS_MT ( GATK4_MERGEBAMALIGNMENT_MT.out.bam )
        ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS_MT.out.versions.first())

        // Marks duplicates
        PICARD_MARKDUPLICATES_MT (PICARD_ADDORREPLACEREADGROUPS_MT.out.bam, genome_fasta, genome_fai )
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES_MT.out.versions.first())

        // Sort bam file
        SAMTOOLS_SORT_MT (PICARD_MARKDUPLICATES_MT.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_SORT_MT.out.versions.first())

        // Index bam file
        SAMTOOLS_INDEX_MT(SAMTOOLS_SORT_MT.out.bam)
        ch_sort_index_bam = SAMTOOLS_SORT_MT.out.bam.join(SAMTOOLS_INDEX_MT.out.bai, by: [0])
        ch_sort_index_bam_intervals_mt = ch_sort_index_bam.combine(intervals_mt)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MT.out.versions.first())

        // Calls variants with Mutect2
        GATK4_MUTECT2_MT (ch_sort_index_bam_intervals_mt, genome_fasta, genome_fai, genome_dict, [], [], [],[])
        ch_versions = ch_versions.mix(GATK4_MUTECT2_MT.out.versions.first())

        // Haplocheck
        // TODO: probably it will be outside this subworkflow as we want to run
        // with the VCF with the variants from the shifted alignment (to solve the mt circularity issue)
        HAPLOCHECK_MT ( GATK4_MUTECT2_MT.out.vcf )
        ch_versions = ch_versions.mix(HAPLOCHECK_MT.out.versions.first())

        // Filter Mutect2 calls
        ch_mutect_vcf = GATK4_MUTECT2_MT.out.vcf.join(GATK4_MUTECT2_MT.out.tbi, by: [0])
        ch_mutect_out = ch_mutect_vcf.join(GATK4_MUTECT2_MT.out.stats, by: [0])
        ch_to_filt = ch_mutect_out.map {
            meta, vcf, tbi, stats ->
                return [meta, vcf, tbi, stats, [], [], [], []]}
        GATK4_FILTERMUTECTCALLS_MT( ch_to_filt,
            genome_fasta,
            genome_fai,
            genome_dict )
        ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS_MT.out.versions.first())

        // Replace within the vcf sample as a sample name with meta.id
        PICARD_RENAMESAMPLEINVCF_MT(GATK4_FILTERMUTECTCALLS_MT.out.vcf)
        ch_versions = ch_versions.mix(PICARD_RENAMESAMPLEINVCF_MT.out.versions.first())

        TABIX_TABIX_MT(PICARD_RENAMESAMPLEINVCF_MT.out.vcf)
        ch_versions = ch_versions.mix(TABIX_TABIX_MT.out.versions.first())


    emit:
        vcf       = PICARD_RENAMESAMPLEINVCF_MT.out.vcf
        tbi       = TABIX_TABIX_MT.out.tbi
        stats     = GATK4_MUTECT2_MT.out.stats
        filt_sats = GATK4_FILTERMUTECTCALLS_MT.out.stats
        txt       = HAPLOCHECK_MT.out.txt
        html      = HAPLOCHECK_MT.out.html
        versions  = ch_versions
}
