//
// A subworkflow to call SNVs by sentieon dnascope with a machine learning model.
//

include { SENTIEON_DNASCOPE             } from '../../../modules/local/sentieon/dnascope'
include { SENTIEON_DNAMODELAPPLY        } from '../../../modules/local/sentieon/dnamodelapply'
include { BCFTOOLS_FILTER               } from '../../../modules/nf-core/bcftools/filter/main'

workflow CALL_SNV_SENTIEON {
	take:
		input         // channel: [ val(meta), bam, bai ]
		fasta         // path: genome.fasta
		fai           // path: genome.fai
		dbsnp         // path: params.known_dbsnp
		dbsnp_index   // path: params.known_dbsnp
        call_interval // path: params.call_interval
        ml_model      // path: params.ml_model

	main:
		ch_versions = Channel.empty()

        SENTIEON_DNASCOPE ( input, fasta, fai, dbsnp, dbsnp_index, call_interval, ml_model )
        ch_vcf      = SENTIEON_DNASCOPE.out.vcf
        ch_index    = SENTIEON_DNASCOPE.out.vcf_index
		ch_versions = ch_versions.mix(SENTIEON_DNASCOPE.out.versions.first())

        if ( ml_model ) {

            ch_vcf_idx = ch_vcf.join( ch_index )

            SENTIEON_DNAMODELAPPLY ( ch_vcf_idx, fasta, fai, ml_model )
            ch_vcf      = SENTIEON_DNAMODELAPPLY.out.vcf
            ch_index    = SENTIEON_DNAMODELAPPLY.out.vcf_index
            ch_versions = ch_versions.mix(SENTIEON_DNAMODELAPPLY.out.versions.first())
        }

        BCFTOOLS_FILTER ( ch_vcf  )
        ch_vcf   = BCFTOOLS_FILTER.out.vcf

	emit:
		vcf		 = ch_vcf
        tabix    = ch_index
        versions = ch_versions.ifEmpty(null)     // channel: [ versions.yml ]
}
