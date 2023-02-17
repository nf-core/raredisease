//
// Align and call MT
//

include { SENTIEON_BWAMEM as SENTIEON_BWAMEM_MT                             } from '../../../modules/local/sentieon/bwamem'
include { BWAMEM2_MEM as BWAMEM2_MEM_MT                                     } from '../../../modules/nf-core/bwamem2/mem/main'
include { GATK4_MERGEBAMALIGNMENT as GATK4_MERGEBAMALIGNMENT_MT             } from '../../../modules/nf-core/gatk4/mergebamalignment/main'
include { PICARD_ADDORREPLACEREADGROUPS as PICARD_ADDORREPLACEREADGROUPS_MT } from '../../../modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_MARKDUPLICATES as PICARD_MARKDUPLICATES_MT                 } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MT                               } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_MT                                 } from '../../../modules/nf-core/samtools/sort/main'
include { HAPLOCHECK as HAPLOCHECK_MT                                       } from '../../../modules/nf-core/haplocheck/main'
include { GATK4_MUTECT2 as GATK4_MUTECT2_MT                                 } from '../../../modules/nf-core/gatk4/mutect2/main'
include { GATK4_FILTERMUTECTCALLS as  GATK4_FILTERMUTECTCALLS_MT            } from '../../../modules/nf-core/gatk4/filtermutectcalls/main'
include { TABIX_TABIX as TABIX_TABIX_MT                                     } from '../../../modules/nf-core/tabix/tabix/main'

workflow ALIGN_AND_CALL_MT {
    take:
        fastq          // channel: [ val(meta), path('*.fastq.gz') ]
        ubam           // channel: [ val(meta), path('*.bam') ]
        index_bwa      // channel: [ /path/to/bwamem2/index/ ]
        index_bwamem2  // channel: [ /path/to/bwamem2/index/ ]
        fasta          // channel: [ genome.fasta ]
        dict           // channel: [ genome.dict ]
        fai            // channel: [ genome.fai ]
        intervals_mt   // channel: [ file(non_control_region.chrM.interval_list) ]

    main:
        ch_versions = Channel.empty()

        BWAMEM2_MEM_MT (fastq , index_bwamem2, true)

        SENTIEON_BWAMEM_MT ( fastq, fasta, fai, index_bwa )

        ch_mt_bam      = Channel.empty().mix(BWAMEM2_MEM_MT.out.bam, SENTIEON_BWAMEM_MT.out.bam)
        ch_fastq_ubam  = ch_mt_bam.join(ubam, by: [0])

        GATK4_MERGEBAMALIGNMENT_MT (ch_fastq_ubam, fasta, dict)

        PICARD_ADDORREPLACEREADGROUPS_MT (GATK4_MERGEBAMALIGNMENT_MT.out.bam)

        PICARD_MARKDUPLICATES_MT (PICARD_ADDORREPLACEREADGROUPS_MT.out.bam, fasta, fai)

        SAMTOOLS_SORT_MT (PICARD_MARKDUPLICATES_MT.out.bam)

        SAMTOOLS_INDEX_MT(SAMTOOLS_SORT_MT.out.bam)
        ch_sort_index_bam        = SAMTOOLS_SORT_MT.out.bam.join(SAMTOOLS_INDEX_MT.out.bai, by: [0])
        ch_sort_index_bam_int_mt = ch_sort_index_bam.combine(intervals_mt)

        GATK4_MUTECT2_MT (ch_sort_index_bam_int_mt, fasta, fai, dict, [], [], [],[])

        // Haplocheck
        // TODO: probably it will be outside this subworkflow as we want to run
        // with the VCF with the variants from the shifted alignment (to solve the mt circularity issue)
        HAPLOCHECK_MT (GATK4_MUTECT2_MT.out.vcf)

        // Filter Mutect2 calls
        ch_mutect_vcf = GATK4_MUTECT2_MT.out.vcf.join(GATK4_MUTECT2_MT.out.tbi, by: [0])
        ch_mutect_out = ch_mutect_vcf.join(GATK4_MUTECT2_MT.out.stats, by: [0])
        ch_to_filt    = ch_mutect_out.map {
                            meta, vcf, tbi, stats ->
                            return [meta, vcf, tbi, stats, [], [], [], []]
                        }

        GATK4_FILTERMUTECTCALLS_MT (ch_to_filt, fasta, fai, dict)

        ch_versions = ch_versions.mix(BWAMEM2_MEM_MT.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_MERGEBAMALIGNMENT_MT.out.versions.first())
        ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS_MT.out.versions.first())
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES_MT.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_SORT_MT.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MT.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_MUTECT2_MT.out.versions.first())
        ch_versions = ch_versions.mix(HAPLOCHECK_MT.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS_MT.out.versions.first())

    emit:
        vcf       = GATK4_FILTERMUTECTCALLS_MT.out.vcf
        tbi       = GATK4_FILTERMUTECTCALLS_MT.out.tbi
        stats     = GATK4_MUTECT2_MT.out.stats
        filt_sats = GATK4_FILTERMUTECTCALLS_MT.out.stats
        txt       = HAPLOCHECK_MT.out.txt
        html      = HAPLOCHECK_MT.out.html
        versions  = ch_versions
}
