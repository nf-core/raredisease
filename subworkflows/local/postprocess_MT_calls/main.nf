//
// Merge and normalize MT variants
//

include { GATK4_MERGEVCFS as GATK4_MERGEVCFS_LIFT_UNLIFT_MT     } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_VARIANTFILTRATION as GATK4_VARIANTFILTRATION_MT } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_MT               } from '../../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX_TABIX_MT                         } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES_MT                 } from '../../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX_TABIX_MT2                        } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_MT                   } from '../../../modules/nf-core/bcftools/merge/main'
include { TABIX_TABIX as TABIX_TABIX_MERGE                      } from '../../../modules/nf-core/tabix/tabix/main'
include { PICARD_LIFTOVERVCF                                    } from '../../../modules/nf-core/picard/liftovervcf/main'
include { BCFTOOLS_ANNOTATE                                     } from '../../../modules/nf-core/bcftools/annotate/main'
include { ADD_VARCALLER_TO_BED                                  } from '../../../modules/local/add_varcallername_to_bed'
include { TABIX_TABIX as TABIX_ANNOTATE                         } from '../../../modules/nf-core/tabix/tabix/main'

workflow POSTPROCESS_MT_CALLS {
    take:
        ch_mt_vcf              // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_mtshift_vcf         // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_genome_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_dictionary   // channel: [mandatory] [ val(meta), path(dict) ]
        ch_genome_fai          // channel: [mandatory] [ val(meta), path(fai) ]
        ch_mtshift_backchain   // channel: [mandatory] [ val(meta), path(backchain) ]
        ch_case_info           // channel: [mandatory] [ val(case_info) ]
        ch_foundin_header      // channel: [mandatory] [ path(header) ]
        ch_genome_chrsizes // channel: [mandatory] [ path(chrsizes) ]

    main:
        ch_versions = Channel.empty()

        // LIFTOVER SHIFTED VCF TO REFERENCE MT POSITIONS
        PICARD_LIFTOVERVCF (
            ch_mtshift_vcf,
            ch_genome_dictionary,
            ch_genome_fasta,
            ch_mtshift_backchain,
        )

        ch_vcfs = ch_mt_vcf
            .join(PICARD_LIFTOVERVCF.out.vcf_lifted, remainder: true)
            .map{ meta, vcf1, vcf2 ->
                [meta, [vcf1, vcf2]]
            }

        GATK4_MERGEVCFS_LIFT_UNLIFT_MT( ch_vcfs, ch_genome_dictionary)

        // Filtering Variants
        GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.vcf
            .join(GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.tbi, failOnMismatch:true, failOnDuplicate:true)
            .set { ch_filt_vcf }
        GATK4_VARIANTFILTRATION_MT (ch_filt_vcf, ch_genome_fasta, ch_genome_fai, ch_genome_dictionary)

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
            .map{ it -> it[1]}
            .toSortedList{a, b -> a.name <=> b.name}
            .toList()
            .set { file_list_vcf }

        TABIX_TABIX_MT2.out.tbi
            .map{ it -> it[1]}
            .toSortedList{a, b -> a.name <=> b.name}
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
            [[:],[]]
        )

        BCFTOOLS_MERGE_MT.out.vcf
            .mix(ch_case_vcf.single)
            .set { ch_addfoundintag_in }

        TABIX_TABIX_MERGE(ch_addfoundintag_in)

        ch_genome_chrsizes.flatten().map{chromsizes ->
            return [[id:'mutect2'], chromsizes]
            }
            .set { ch_varcallerinfo }

        ADD_VARCALLER_TO_BED (ch_varcallerinfo).gz_tbi
            .map{meta,bed,tbi -> return [bed, tbi]}
            .set{ch_varcallerbed}

        ch_addfoundintag_in
            .join(TABIX_TABIX_MERGE.out.tbi)
            .combine(ch_varcallerbed)
            .set { ch_annotate_in }

        BCFTOOLS_ANNOTATE(ch_annotate_in,ch_foundin_header)

        TABIX_ANNOTATE(BCFTOOLS_ANNOTATE.out.vcf)

        ch_versions = ch_versions.mix(PICARD_LIFTOVERVCF.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION_MT.out.versions.first())
        ch_versions = ch_versions.mix(SPLIT_MULTIALLELICS_MT.out.versions.first())
        ch_versions = ch_versions.mix(REMOVE_DUPLICATES_MT.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_MERGE_MT.out.versions)
        ch_versions = ch_versions.mix(ADD_VARCALLER_TO_BED.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions)
        ch_versions = ch_versions.mix(TABIX_ANNOTATE.out.versions)
        ch_versions = ch_versions.mix(TABIX_TABIX_MERGE.out.versions)
        ch_versions = ch_versions.mix(TABIX_TABIX_MT.out.versions)
        ch_versions = ch_versions.mix(TABIX_TABIX_MT2.out.versions)

    emit:
        vcf       = BCFTOOLS_ANNOTATE.out.vcf   // channel: [ val(meta), path(vcf) ]
        tbi       = TABIX_ANNOTATE.out.tbi      // channel: [ val(meta), path(tbi) ]
        versions  = ch_versions                 // channel: [ path(versions.yml) ]
}
