//
// A subworkflow to call SNVs by sentieon dnascope with a machine learning model.
//

include { SENTIEON_DNASCOPE                          } from '../../../modules/nf-core/sentieon/dnascope/main'
include { SENTIEON_DNAMODELAPPLY                     } from '../../../modules/nf-core/sentieon/dnamodelapply/main'
include { BCFTOOLS_MERGE                             } from '../../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_SEN   } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES_SEN     } from '../../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX_SEN                   } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_BCFTOOLS              } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_FILTER as BCF_FILTER_ONE          } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as BCF_FILTER_TWO          } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_ANNOTATE                          } from '../../../modules/nf-core/bcftools/annotate/main'
include { TABIX_TABIX as TABIX_ANNOTATE              } from '../../../modules/nf-core/tabix/tabix/main'
include { ADD_VARCALLER_TO_BED                       } from '../../../modules/local/add_varcallername_to_bed'

workflow CALL_SNV_SENTIEON {
    take:
        ch_bam_bai         // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_fasta    // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai      // channel: [mandatory] [ val(meta), path(fai) ]
        ch_dbsnp           // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_dbsnp_index     // channel: [mandatory] [ val(meta), path(tbi) ]
        ch_call_interval   // channel: [mandatory] [ val(meta), path(interval) ]
        ch_ml_model        // channel: [mandatory] [ val(meta), path(model) ]
        ch_case_info       // channel: [mandatory] [ val(case_info) ]
        ch_pcr_indel_model // channel: [optional] [ val(sentieon_dnascope_pcr_indel_model) ]
        ch_foundin_header  // channel: [mandatory] [ path(header) ]
        ch_genome_chrsizes // channel: [mandatory] [ path(chrsizes) ]

    main:
        ch_versions = Channel.empty()

        // Combine bam and intervals
        bam_bai_intervals = ch_bam_bai.combine(ch_call_interval)
            .map{
                meta, bam, bai, meta2, interval -> [meta, bam, bai, interval]
            }

        SENTIEON_DNASCOPE(
            bam_bai_intervals,
            ch_genome_fasta,
            ch_genome_fai,
            ch_dbsnp,
            ch_dbsnp_index,
            ch_ml_model,
            ch_pcr_indel_model,
            'VARIANT',
            false)

        ch_dnamodelapply_in = SENTIEON_DNASCOPE.out.vcf.join(SENTIEON_DNASCOPE.out.vcf_tbi)

        SENTIEON_DNAMODELAPPLY ( ch_dnamodelapply_in, ch_genome_fasta, ch_genome_fai, ch_ml_model )

        BCF_FILTER_ONE (SENTIEON_DNAMODELAPPLY.out.vcf )

        BCF_FILTER_TWO ( BCF_FILTER_ONE.out.vcf )

        TABIX_BCFTOOLS ( BCF_FILTER_TWO.out.vcf )

        BCF_FILTER_TWO.out.vcf.join(TABIX_BCFTOOLS.out.tbi, failOnMismatch:true, failOnDuplicate:true)
            .map { meta,vcf,tbi -> return [vcf,tbi] }
            .set { ch_vcf_idx }

        ch_case_info
            .combine(ch_vcf_idx)
            .groupTuple()
            .branch{                                                                                                    // branch the channel into multiple channels (single, multiple) depending on size of list
                single: it[1].size() == 1
                multiple: it[1].size() > 1
            }
            .set{ ch_vcf_idx_merge_in }

        BCFTOOLS_MERGE(ch_vcf_idx_merge_in.multiple, ch_genome_fasta, ch_genome_fai, [])

        ch_split_multi_in = BCFTOOLS_MERGE.out.merged_variants
                    .map{meta, bcf ->
                        return [meta, bcf, []]}

        ch_vcf_idx_case =  ch_vcf_idx_merge_in.single.mix(ch_split_multi_in)

        SPLIT_MULTIALLELICS_SEN(ch_vcf_idx_case, ch_genome_fasta)

        ch_remove_dup_in = SPLIT_MULTIALLELICS_SEN.out.vcf
                            .map{meta, vcf ->
                                    return [meta, vcf, []]}

        REMOVE_DUPLICATES_SEN(ch_remove_dup_in, ch_genome_fasta)

        TABIX_SEN(REMOVE_DUPLICATES_SEN.out.vcf)

        ch_genome_chrsizes.flatten().map{chromsizes ->
            return [[id:'sentieon_dnascope'], chromsizes]
            }
            .set { ch_varcallerinfo }

        ADD_VARCALLER_TO_BED (ch_varcallerinfo).gz_tbi
            .map{meta,bed,tbi -> return [bed, tbi]}
            .set{ch_varcallerbed}

        REMOVE_DUPLICATES_SEN.out.vcf
            .join(TABIX_SEN.out.tbi)
            .combine(ch_varcallerbed)
            .combine(ch_foundin_header)
            .set { ch_annotate_in }

        BCFTOOLS_ANNOTATE(ch_annotate_in)

        TABIX_ANNOTATE(BCFTOOLS_ANNOTATE.out.vcf)

        ch_versions = ch_versions.mix(SENTIEON_DNASCOPE.out.versions.first())
        ch_versions = ch_versions.mix(SENTIEON_DNAMODELAPPLY.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions.first())
        ch_versions = ch_versions.mix(SPLIT_MULTIALLELICS_SEN.out.versions.first())
        ch_versions = ch_versions.mix(REMOVE_DUPLICATES_SEN.out.versions.first())
        ch_versions = ch_versions.mix(TABIX_SEN.out.versions.first())
        ch_versions = ch_versions.mix(BCF_FILTER_ONE.out.versions.first())
        ch_versions = ch_versions.mix(ADD_VARCALLER_TO_BED.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions)
        ch_versions = ch_versions.mix(TABIX_ANNOTATE.out.versions)

    emit:
        vcf      = BCFTOOLS_ANNOTATE.out.vcf  // channel: [ val(meta), path(vcf) ]
        tabix    = TABIX_ANNOTATE.out.tbi     // channel: [ val(meta), path(tbi) ]
        versions = ch_versions                // channel: [ path(versions.yml) ]
}
