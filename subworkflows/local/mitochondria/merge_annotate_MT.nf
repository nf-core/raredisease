//
// Merge and annotate MT
//

include { GATK4_MERGEVCFS as GATK4_MERGEVCFS_LIFT_UNLIFT_MT     } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_VARIANTFILTRATION as GATK4_VARIANTFILTRATION_MT } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_MT               } from '../../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX_TABIX_MT                         } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES_MT                 } from '../../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX_TABIX_MT2                        } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_MT                   } from '../../../modules/nf-core/bcftools/merge/main'
include { TABIX_TABIX as TABIX_TABIX_MT3                        } from '../../../modules/nf-core/tabix/tabix/main'
include { ENSEMBLVEP as ENSEMBLVEP_MT                           } from '../../../modules/local/ensemblvep/main'
include { HAPLOGREP2_CLASSIFY as HAPLOGREP2_CLASSIFY_MT         } from '../../../modules/nf-core/haplogrep2/classify/main'
include { VCFANNO as VCFANNO_MT                                 } from '../../../modules/nf-core/vcfanno/main'
include { TABIX_BGZIPTABIX as ZIP_TABIX_HMTNOTE                 } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { HMTNOTE_ANNOTATE as HMTNOTE_ANNOTATE                  } from '../../../modules/nf-core/hmtnote/annotate/main'

workflow MERGE_ANNOTATE_MT {
    take:
        ch_vcf1                // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_vcf2                // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_genome_fasta        // channel: [mandatory] [ path(fasta) ]
        ch_genome_dict_meta    // channel: [mandatory] [ val(meta), path(dict) ]
        ch_genome_dict_no_meta // channel: [mandatory] [ path(dict) ]
        ch_genome_fai          // channel: [mandatory] [ path(fai) ]
        ch_vcfanno_resources   // channel: [mandatory] [ path(resources) ]
        ch_vcfanno_toml        // channel: [mandatory] [ path(toml) ]
        val_vep_genome         // string:  [mandatory] GRCh37 or GRCh38
        val_vep_cache_version  // string:  [mandatory] 107
        ch_vep_cache           // channel: [mandatory] [ path(cache) ]
        ch_case_info           // channel: [mandatory] [ val(case_info) ]

    main:
        ch_versions = Channel.empty()

        ch_vcfs = ch_vcf1
            .join(ch_vcf2, remainder: true)
            .map{ meta, vcf1, vcf2 ->
            [meta, [vcf1, vcf2]]
        }
        GATK4_MERGEVCFS_LIFT_UNLIFT_MT( ch_vcfs, ch_genome_dict_meta)

        // Filtering Variants
        GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.vcf
            .join(GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.tbi, failOnMismatch:true, failOnDuplicate:true)
            .set { ch_filt_vcf }
        GATK4_VARIANTFILTRATION_MT (ch_filt_vcf, ch_genome_fasta, ch_genome_fai, ch_genome_dict_no_meta)

        // Spliting multiallelic calls
        GATK4_VARIANTFILTRATION_MT.out.vcf
            .join(GATK4_VARIANTFILTRATION_MT.out.tbi, failOnMismatch:true, failOnDuplicate:true)
            .set { ch_in_split }
        SPLIT_MULTIALLELICS_MT (ch_in_split, ch_genome_fasta)
        TABIX_TABIX_MT(SPLIT_MULTIALLELICS_MT.out.vcf)

        // Removing duplicates and merging if there is more than one sample
        SPLIT_MULTIALLELICS_MT.out.vcf
            .join(TABIX_TABIX_MT.out.tbi, failOnMismatch:true, failOnDuplicate:true)
            .set { ch_in_remdup }
        REMOVE_DUPLICATES_MT(ch_in_remdup, ch_genome_fasta)
        TABIX_TABIX_MT2(REMOVE_DUPLICATES_MT.out.vcf)

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

        ch_case_info
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
            ch_genome_fasta,
            ch_genome_fai)
        ch_merged_vcf = BCFTOOLS_MERGE_MT.out.merged_variants

        ch_in_vep = ch_merged_vcf.mix(ch_case_vcf.single)

        // Annotating with ensembl Vep
        ENSEMBLVEP_MT( ch_in_vep,
            val_vep_genome,
            "homo_sapiens",
            val_vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta,
            [])

        // Running vcfanno
        TABIX_TABIX_MT3(ENSEMBLVEP_MT.out.vcf_gz)
        ch_in_vcfanno = ENSEMBLVEP_MT.out.vcf_gz.join(TABIX_TABIX_MT3.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        VCFANNO_MT(ch_in_vcfanno, ch_vcfanno_toml, [], ch_vcfanno_resources)

        // HMTNOTE ANNOTATE
        HMTNOTE_ANNOTATE(VCFANNO_MT.out.vcf)
        ZIP_TABIX_HMTNOTE(HMTNOTE_ANNOTATE.out.vcf)

        // Prepare output
        ch_vcf_out = ZIP_TABIX_HMTNOTE.out.gz_tbi.map{meta, vcf, tbi -> return [meta, vcf] }
        ch_tbi_out = ZIP_TABIX_HMTNOTE.out.gz_tbi.map{meta, vcf, tbi -> return [meta, tbi] }

        // Running haplogrep2
        HAPLOGREP2_CLASSIFY_MT(ch_in_vep, "vcf.gz")

        ch_versions = ch_versions.mix(GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION_MT.out.versions.first())
        ch_versions = ch_versions.mix(SPLIT_MULTIALLELICS_MT.out.versions.first())
        ch_versions = ch_versions.mix(REMOVE_DUPLICATES_MT.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_MERGE_MT.out.versions)
        ch_versions = ch_versions.mix(ENSEMBLVEP_MT.out.versions)
        ch_versions = ch_versions.mix(VCFANNO_MT.out.versions)
        ch_versions = ch_versions.mix(HMTNOTE_ANNOTATE.out.versions)
        ch_versions = ch_versions.mix(HAPLOGREP2_CLASSIFY_MT.out.versions)

    emit:
        haplog    = HAPLOGREP2_CLASSIFY_MT.out.txt // channel: [ val(meta), path(txt) ]
        vcf       = ch_vcf_out                     // channel: [ val(meta), path(vcf) ]
        tbi       = ch_tbi_out                     // channel: [ val(meta), path(tbi) ]
        report    = ENSEMBLVEP_MT.out.report       // channel: [ path(html) ]
        versions  = ch_versions                    // channel: [ path(versions.yml) ]
}
