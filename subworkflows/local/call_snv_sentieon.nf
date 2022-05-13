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
        	ch_dnascope_vcf = Channel.empty()
        	ch_dnamodelapply_vcf = Channel.empty()

	        SENTIEON_DNASCOPE ( input, fasta, fai, dbsnp, dbsnp_index, ml_model )
        	ch_versions = ch_versions.mix(SENTIEON_DNASCOPE.out.versions)

        	SENTIEON_DNASCOPE.out
            	.vcf
            	.set { ch_dnascope_vcf }

        	SENTIEON_DNAMODELAPPLY ( ch_dnascope_vcf, fasta, fai, ml_model )
        	ch_versions = ch_versions.mix(SENTIEON_DNAMODELAPPLY.out.versions)
        	ch_dnamodelapply_vcf = SENTIEON_DNAMODELAPPLY.out.vcf

	emit:
		dnascope_vcf                = ch_dnascope_vcf.ifEmpty(null)
        	dnamodelapply_vcf           = ch_dnamodelapply_vcf.ifEmpty(null)
        	versions                    = ch_versions.ifEmpty(null)             // channel: [ versions.yml ]
}
