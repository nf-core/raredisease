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
include { TABIX_TABIX as TABIX_TABIX_MERGE                      } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TABIX_MT3                        } from '../../../modules/nf-core/tabix/tabix/main'
include { ENSEMBLVEP as ENSEMBLVEP_MT                           } from '../../../modules/local/ensemblvep/main'
include { HAPLOGREP2_CLASSIFY as HAPLOGREP2_CLASSIFY_MT         } from '../../../modules/nf-core/haplogrep2/classify/main'
include { VCFANNO as VCFANNO_MT                                 } from '../../../modules/nf-core/vcfanno/main'
include { ANNOTATE_CADD                                         } from '../annotation/annotate_cadd'
include { TABIX_BGZIPTABIX as ZIP_TABIX_HMTNOTE                 } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { HMTNOTE_ANNOTATE                                      } from '../../../modules/nf-core/hmtnote/annotate/main'

workflow MERGE_ANNOTATE_MT {
    take:
        ch_vcf1                // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_vcf2                // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_cadd_header         // channel: [mandatory] [ path(txt) ]
        ch_cadd_resources      // channel: [mandatory] [ path(annotation) ]
        ch_genome_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_dict         // channel: [mandatory] [ val(meta), path(dict) ]
        ch_genome_fai          // channel: [mandatory] [ val(meta), path(fai) ]
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
        GATK4_MERGEVCFS_LIFT_UNLIFT_MT( ch_vcfs, ch_genome_dict)

        // Filtering Variants
        GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.vcf
            .join(GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.tbi, failOnMismatch:true, failOnDuplicate:true)
            .set { ch_filt_vcf }
        GATK4_VARIANTFILTRATION_MT (ch_filt_vcf, ch_genome_fasta, ch_genome_fai, ch_genome_dict)

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
            ch_genome_fasta,
            ch_genome_fai,
            []
        )

        BCFTOOLS_MERGE_MT.out.merged_variants
            .mix(ch_case_vcf.single)
            .set { ch_annotation_in }

        TABIX_TABIX_MERGE(ch_annotation_in)

        // Annotating with CADD
        ANNOTATE_CADD (
            ch_annotation_in,
            TABIX_TABIX_MERGE.out.tbi,
            ch_cadd_header,
            ch_cadd_resources
        )

        // Pick input for vep
        ch_annotation_in
            .combine(ANNOTATE_CADD.out.vcf.ifEmpty("null"))
            .branch { it  ->
                merged: it[2].equals("null")
                    return [it[0], it[1]]
                cadd: !(it[2].equals("null"))
                    return [it[2], it[3]]
            }
            .set { ch_for_mix }
        ch_vep_in = ch_for_mix.merged.mix(ch_for_mix.cadd)

        // Annotating with ensembl Vep
        ENSEMBLVEP_MT(
            ch_vep_in,
            ch_genome_fasta,
            val_vep_genome,
            "homo_sapiens",
            val_vep_cache_version,
            ch_vep_cache,
            []
        )

        // Running vcfanno
        TABIX_TABIX_MT3(ENSEMBLVEP_MT.out.vcf_gz)
        ENSEMBLVEP_MT.out.vcf_gz
            .join(TABIX_TABIX_MT3.out.tbi, failOnMismatch:true, failOnDuplicate:true)
            .map { meta, vcf, tbi -> return [meta, vcf, tbi, []]}
            .set { ch_in_vcfanno }

        VCFANNO_MT(ch_in_vcfanno, ch_vcfanno_toml, [], ch_vcfanno_resources)

        // HMTNOTE ANNOTATE
        HMTNOTE_ANNOTATE(VCFANNO_MT.out.vcf)
        HMTNOTE_ANNOTATE.out.vcf.map{meta, vcf ->
            return [meta, WorkflowRaredisease.replaceSpacesInInfoColumn(vcf, vcf.parent.toString(), vcf.baseName)]
            }
            .set { ch_hmtnote_reformatted }
        ZIP_TABIX_HMTNOTE(ch_hmtnote_reformatted)

        // Prepare output
        ch_vcf_out = ZIP_TABIX_HMTNOTE.out.gz_tbi.map{meta, vcf, tbi -> return [meta, vcf] }
        ch_tbi_out = ZIP_TABIX_HMTNOTE.out.gz_tbi.map{meta, vcf, tbi -> return [meta, tbi] }

        // Running haplogrep2
        HAPLOGREP2_CLASSIFY_MT(ch_vep_in, "vcf.gz")

        ch_versions = ch_versions.mix(GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION_MT.out.versions.first())
        ch_versions = ch_versions.mix(SPLIT_MULTIALLELICS_MT.out.versions.first())
        ch_versions = ch_versions.mix(REMOVE_DUPLICATES_MT.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_MERGE_MT.out.versions)
        ch_versions = ch_versions.mix(ANNOTATE_CADD.out.versions)
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
