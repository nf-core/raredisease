//
// A subworkflow to call SNVs by sentieon dnascope with a machine learning model.
//

include { SENTIEON_DNASCOPE        }   from '../../modules/local/sentieon/dnascope'
include { SENTIEON_DNAMODELAPPLY   }   from '../../modules/local/sentieon/dnamodelapply'

workflow CALL_SNV_SENTIEON {
	take:
		input           // channel: [ val(meta), bam, bai ]
		fasta           // path: genome.fasta
		fai             // path: genome.fai
		known_dbsnp     // path: params.known_dbsnp
		known_dbsnp_tbi // path: params.known_dbsnp
        ml_model        // path: params.ml_model

	main:
		ch_versions = Channel.empty()

        SENTIEON_DNASCOPE ( input, fasta, fai, dbsnp, dbsnp_index, ml_model )
        ch_vcf = SENTIEON_DNASCOPE.out.vcf
        ch_vcf_index = SENTIEON_DNASCOPE.out.vcf_index
		ch_versions = ch_versions.mix(SENTIEON_DNASCOPE.out.versions)

        if ( ml_model ) {

            ch_vcf.
            .join( ch_vcf_index )
            .set { ch_vcf_idx }

            SENTIEON_DNAMODELAPPLY ( ch_vcf_idx, fasta, fai, ml_model )
            ch_vcf = SENTIEON_DNAMODELAPPLY.out.vcf
            ch_vcf_index = SENTIEON_DNAMODELAPPLY.out.vcf_index
            ch_versions = ch_versions.mix(SENTIEON_DNAMODELAPPLY.out.versions)
        }

	emit:
		vcf			= ch_vcf
        vcf_index   = ch_vcf_index
        versions	= ch_versions.ifEmpty(null)     // channel: [ versions.yml ]
}
