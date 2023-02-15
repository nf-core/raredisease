//
// A subworkflow to call SNVs by sentieon dnascope with a machine learning model.
//

include { SENTIEON_DNASCOPE                        } from '../../../modules/local/sentieon/dnascope'
include { SENTIEON_DNAMODELAPPLY                   } from '../../../modules/local/sentieon/dnamodelapply'
include { BCFTOOLS_MERGE                           } from '../../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_SEN } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES_SEN   } from '../../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX_SEN                 } from '../../../modules/nf-core/tabix/tabix/main'

workflow CALL_SNV_SENTIEON {
	take:
		input         // channel: [ val(meta), bam, bai ]
		fasta         // path: genome.fasta
		fai           // path: genome.fai
		dbsnp         // path: params.known_dbsnp
		dbsnp_index   // path: params.known_dbsnp
        call_interval // path: params.call_interval
        ml_model      // path: params.ml_model
        case_info     // channel: [ case_id ]

	main:
		ch_versions = Channel.empty()

        SENTIEON_DNASCOPE ( input, fasta, fai, dbsnp, dbsnp_index, call_interval, ml_model )
        ch_vcf      = SENTIEON_DNASCOPE.out.vcf
        ch_index    = SENTIEON_DNASCOPE.out.vcf_index

        ch_vcf_idx  = ch_vcf.join( ch_index )

        SENTIEON_DNAMODELAPPLY ( ch_vcf_idx, fasta, fai, ml_model )

        if (params.ml_model) {
            ch_vcf      = SENTIEON_DNAMODELAPPLY.out.vcf
            ch_index    = SENTIEON_DNAMODELAPPLY.out.vcf_index
        }

        ch_vcf.join(ch_index)
            .map { meta,vcf,tbi -> return [vcf,tbi] }
            .set { ch_vcf_idx }

        case_info
            .combine(ch_vcf_idx)
            .groupTuple()
            .branch{                                                                                                    // branch the channel into multiple channels (single, multiple) depending on size of list
                single: it[1].size() == 1
                multiple: it[1].size() > 1
            }
            .set{ ch_vcf_idx_merge_in }

        BCFTOOLS_MERGE(ch_vcf_idx_merge_in.multiple,[],fasta,fai)
        ch_split_multi_in = BCFTOOLS_MERGE.out.merged_variants
                    .map{meta, bcf ->
                        return [meta, bcf, []]}
        ch_vcf_idx_case =  ch_vcf_idx_merge_in.single.mix(ch_split_multi_in)
        SPLIT_MULTIALLELICS_SEN(ch_vcf_idx_case, fasta)

        ch_remove_dup_in = SPLIT_MULTIALLELICS_SEN.out.vcf
                            .map{meta, vcf ->
                                    return [meta, vcf, []]}
        REMOVE_DUPLICATES_SEN(ch_remove_dup_in, fasta)
        TABIX_SEN(REMOVE_DUPLICATES_SEN.out.vcf)

		ch_versions = ch_versions.mix(SENTIEON_DNASCOPE.out.versions.first())
        ch_versions = ch_versions.mix(SENTIEON_DNAMODELAPPLY.out.versions.first())
		ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions.first())
		ch_versions = ch_versions.mix(SPLIT_MULTIALLELICS_SEN.out.versions.first())
		ch_versions = ch_versions.mix(REMOVE_DUPLICATES_SEN.out.versions.first())
		ch_versions = ch_versions.mix(TABIX_SEN.out.versions.first())

	emit:
		vcf		 = REMOVE_DUPLICATES_SEN.out.vcf
        tabix    = TABIX_SEN.out.tbi
        versions = ch_versions.ifEmpty(null)     // channel: [ versions.yml ]
}
