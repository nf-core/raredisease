//
// Merge and annotate MT
//

include { GATK4_MERGEVCFS as GATK4_MERGEVCFS_LIFT_UNLIFT_MT      } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_VARIANTFILTRATION as GATK4_VARIANTFILTRATION_MT  } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_MT                } from '../../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX_TABIX_MT                          } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES_MT                  } from '../../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX_TABIX_MT2                         } from '../../../modules/nf-core/tabix/tabix/main'
include { CHANGE_NAME as CHANGE_NAME_VCF_MT                      } from '../../../modules/local/change_name'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_MT                    } from '../../../modules/nf-core/bcftools/merge/main'
include { HMTNOTE as HMTNOTE_MT                                  } from '../../../modules/nf-core/hmtnote/main'
include { TABIX_TABIX as TABIX_TABIX_MT3                         } from '../../../modules/nf-core/tabix/tabix/main'
include { VCFANNO as VCFANNO_MT                                  } from '../../../modules/nf-core/vcfanno/main'
include { TABIX_TABIX as TABIX_TABIX_MT4                         } from '../../../modules/nf-core/tabix/tabix/main'
include { ENSEMBLVEP as ENSEMBLVEP_MT                            } from '../../../modules/local/ensemblvep/main'
include { HAPLOGREP2_CLASSIFY as HAPLOGREP2_CLASSIFY_MT          } from '../../../modules/nf-core/haplogrep2/classify/main'

workflow MERGE_ANNOTATE_MT {
    take:
        vcf1                // channel: [ val(meta), path('*.vcf.gz') ]
        vcf2                // channel: [ val(meta), path('*.vcf.gz') ]
        genome_fasta        // channel: [ genome.fasta ]
        genome_dict_meta    // channel: [ genome.dict ]
        genome_dict_no_meta // channel: [ genome.dict ]
        genome_fai          // channel: [ genome.fai ]
        vcfanno_resources
        vcfanno_lua
        vcfanno_toml
        vep_genome
        vep_cache_version
        vep_cache
        case_info           // channel: [ val(case_info) ]

    main:
        ch_versions = Channel.empty()

        ch_vcfs = vcf1
            .join(vcf2, remainder: true)
            .map{ meta, vcf1, vcf2 ->
            [meta, [vcf1, vcf2]]
        }

        GATK4_MERGEVCFS_LIFT_UNLIFT_MT( ch_vcfs, genome_dict_meta)
        ch_versions = ch_versions.mix(GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.versions.first())

        // Filtering Variants
        ch_filt_vcf = GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.vcf.join(GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.tbi, by:[0])
        GATK4_VARIANTFILTRATION_MT(ch_filt_vcf,
            genome_fasta,
            genome_fai,
            genome_dict_no_meta )
        ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION_MT.out.versions.first())

        // Spliting multiallelic calls
        ch_in_split = GATK4_VARIANTFILTRATION_MT.out.vcf.join(GATK4_VARIANTFILTRATION_MT.out.tbi, by:[0])
        SPLIT_MULTIALLELICS_MT (ch_in_split, genome_fasta)
        ch_versions = ch_versions.mix(SPLIT_MULTIALLELICS_MT.out.versions.first())

        TABIX_TABIX_MT(SPLIT_MULTIALLELICS_MT.out.vcf)
        ch_in_remdup = SPLIT_MULTIALLELICS_MT.out.vcf.join(TABIX_TABIX_MT.out.tbi)

        // Removing duplicates and merging if there is more than one sample
        REMOVE_DUPLICATES_MT(ch_in_remdup, genome_fasta)
        TABIX_TABIX_MT2(REMOVE_DUPLICATES_MT.out.vcf)
        ch_versions = ch_versions.mix(REMOVE_DUPLICATES_MT.out.versions.first())

        REMOVE_DUPLICATES_MT.out.vcf
            .collect{it[1]}
            .ifEmpty([])
            .toList()
            .set { file_list_vcf }

        TABIX_TABIX_MT2.out.tbi
            .collect{it[1]}
            .ifEmpty([])
            .toList()
            .set { file_list_tbi }

        case_info
            .combine(file_list_vcf)
            .combine(file_list_tbi)
            .set { ch_rem_dup_vcf_tbi }

        ch_rem_dup_vcf_tbi.branch {
            meta, vcf, tbi ->
                single: vcf.size() == 1
                    return [meta, vcf]
                multiple: vcf.size() > 1
                    return [meta, vcf, tbi]
            }.set { ch_case_vcf }

        BCFTOOLS_MERGE_MT( ch_case_vcf.multiple,
            [],
            genome_fasta,
            genome_fai)
        ch_merged_vcf = BCFTOOLS_MERGE_MT.out.merged_variants
        ch_versions = ch_versions.mix(BCFTOOLS_MERGE_MT.out.versions)

        CHANGE_NAME_VCF_MT(ch_case_vcf.single)
        ch_vcf_changed_name = CHANGE_NAME_VCF_MT.out.file
        ch_versions = ch_versions.mix(CHANGE_NAME_VCF_MT.out.versions)

        ch_annot = ch_merged_vcf.mix(ch_vcf_changed_name)

        // Annotating with Hmtnote
        //HMTNOTE_MT(ch_in_vep)
        //ch_versions = ch_versions.mix(HMTNOTE_MT.out.versions.first())

        // Annotating with VCFANNO
        TABIX_TABIX_MT3(ch_annot)
        ch_in_vcfanno=ch_annot.join(TABIX_TABIX_MT3.out.tbi, by:[0])
        VCFANNO_MT(ch_in_vcfanno, vcfanno_toml, vcfanno_lua, vcfanno_resources)
        ch_versions = ch_versions.mix(VCFANNO_MT.out.versions)

        // Annotating with ensembl Vep

        ENSEMBLVEP_MT( VCFANNO_MT.out.vcf,
            vep_genome,
            "homo_sapiens",
            vep_cache_version,
            vep_cache,
            genome_fasta,
            [])
        ch_versions = ch_versions.mix(ENSEMBLVEP_MT.out.versions)

        // Running haplogrep2
        TABIX_TABIX_MT4(ENSEMBLVEP_MT.out.vcf_gz)
        HAPLOGREP2_CLASSIFY_MT(VCFANNO_MT.out.vcf, "vcf.gz")
        ch_versions = ch_versions.mix(HAPLOGREP2_CLASSIFY_MT.out.versions)

    emit:
        haplog   = HAPLOGREP2_CLASSIFY_MT.out.txt
        vcf      = ENSEMBLVEP_MT.out.vcf_gz
        tbi      = TABIX_TABIX_MT4.out.tbi
        report   = ENSEMBLVEP_MT.out.report
        versions = ch_versions // channel: [ versions.yml ]
}
