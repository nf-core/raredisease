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
include { MT_DELETION                                                       } from '../../../modules/local/mt_deletion_script'
include { EKLIPSE as EKLIPSE_MT                                             } from '../../../modules/nf-core/eklipse/main'

workflow ALIGN_AND_CALL_MT {
    take:
        ch_fastq        // channel: [mandatory] [ val(meta), [ path(reads) ] ]
        ch_ubam         // channel: [mandatory] [ val(meta), path(bam) ]
        ch_bwaindex     // channel: [mandatory for sentieon] [ val(meta), path(index) ]
        ch_bwamem2index // channel: [mandatory for bwamem2] [ val(meta), path(index) ]
        ch_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_dict         // channel: [mandatory] [ val(meta), path(dict) ]
        ch_fai          // channel: [mandatory] [ val(meta), path(fai) ]
        ch_intervals    // channel: [mandatory] [ path(interval_list) ]

    main:
        ch_versions = Channel.empty()

        BWAMEM2_MEM_MT (ch_fastq, ch_bwamem2index, true)

        SENTIEON_BWAMEM_MT ( ch_fastq, ch_fasta, ch_fai, ch_bwaindex )

        Channel.empty()
            .mix(BWAMEM2_MEM_MT.out.bam, SENTIEON_BWAMEM_MT.out.bam)
            .join(ch_ubam, failOnMismatch:true, failOnDuplicate:true)
            .set {ch_bam_ubam}

        GATK4_MERGEBAMALIGNMENT_MT (ch_bam_ubam, ch_fasta, ch_dict)

        PICARD_ADDORREPLACEREADGROUPS_MT (GATK4_MERGEBAMALIGNMENT_MT.out.bam)

        PICARD_MARKDUPLICATES_MT (PICARD_ADDORREPLACEREADGROUPS_MT.out.bam, ch_fasta, ch_fai)

        SAMTOOLS_SORT_MT (PICARD_MARKDUPLICATES_MT.out.bam)

        SAMTOOLS_INDEX_MT(SAMTOOLS_SORT_MT.out.bam)
        ch_sort_index_bam        = SAMTOOLS_SORT_MT.out.bam.join(SAMTOOLS_INDEX_MT.out.bai, failOnMismatch:true, failOnDuplicate:true)
        ch_sort_index_bam_int_mt = ch_sort_index_bam.combine(ch_intervals)

        EKLIPSE_MT(ch_sort_index_bam,[])

        MT_DELETION(ch_sort_index_bam, ch_fasta)

        GATK4_MUTECT2_MT (ch_sort_index_bam_int_mt, ch_fasta, ch_fai, ch_dict, [], [], [],[])

        HAPLOCHECK_MT (GATK4_MUTECT2_MT.out.vcf)

        // Filter Mutect2 calls
        ch_mutect_vcf = GATK4_MUTECT2_MT.out.vcf.join(GATK4_MUTECT2_MT.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        ch_mutect_out = ch_mutect_vcf.join(GATK4_MUTECT2_MT.out.stats, failOnMismatch:true, failOnDuplicate:true)
        ch_to_filt    = ch_mutect_out.map {
                            meta, vcf, tbi, stats ->
                            return [meta, vcf, tbi, stats, [], [], [], []]
                        }

        GATK4_FILTERMUTECTCALLS_MT (ch_to_filt, ch_fasta, ch_fai, ch_dict)

        ch_versions = ch_versions.mix(BWAMEM2_MEM_MT.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_MERGEBAMALIGNMENT_MT.out.versions.first())
        ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS_MT.out.versions.first())
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES_MT.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_SORT_MT.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MT.out.versions.first())
        ch_versions = ch_versions.mix(EKLIPSE_MT.out.versions.first())
        ch_versions = ch_versions.mix(MT_DELETION.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_MUTECT2_MT.out.versions.first())
        ch_versions = ch_versions.mix(HAPLOCHECK_MT.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS_MT.out.versions.first())

    emit:
        vcf            = GATK4_FILTERMUTECTCALLS_MT.out.vcf   // channel: [ val(meta), path(vcf) ]
        tbi            = GATK4_FILTERMUTECTCALLS_MT.out.tbi   // channel: [ val(meta), path(tbi) ]
        stats          = GATK4_MUTECT2_MT.out.stats           // channel: [ val(meta), path(stats) ]
        filt_stats     = GATK4_FILTERMUTECTCALLS_MT.out.stats // channel: [ val(meta), path(tsv) ]
        eklipse_del    = EKLIPSE_MT.out.deletions             // channel: [ val(meta), path(csv) ]
        eklipse_genes  = EKLIPSE_MT.out.genes                 // channel: [ val(meta), path(csv) ]
        eklipse_circos = EKLIPSE_MT.out.circos                // channel: [ val(meta), path(png) ]
        txt            = HAPLOCHECK_MT.out.txt                // channel: [ val(meta), path(txt) ]
        html           = HAPLOCHECK_MT.out.html               // channel: [ val(meta), path(html) ]
        mt_del_result  = MT_DELETION.out.mt_del_result        // channel: [ val(meta), path(txt) ]
        versions       = ch_versions                          // channel: [ path(versions.yml) ]
}
