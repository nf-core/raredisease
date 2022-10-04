//
// A subworkflow to call SNVs by sentieon dnascope with a machine learning model.
//

include { ADD_MOST_SEVERE_CSQ as ADD_MOST_SEVERE_CSQ_CLIN } from '../../modules/local/add_most_severe_consequence'
include { ADD_MOST_SEVERE_CSQ as ADD_MOST_SEVERE_CSQ_RES  } from '../../modules/local/add_most_severe_consequence'

workflow ANNOTATE_CSQ {
	take:
		clinical              // channel: [ val(meta), vcf ]
		research              // channel: [ val(meta), vcf ]
		variant_consequences  // path: consequences.txt

	main:
		ch_versions = Channel.empty()

        ADD_MOST_SEVERE_CSQ_CLIN (
            clinical,
            variant_consequences
        )
        ADD_MOST_SEVERE_CSQ_RES (
            research,
            variant_consequences
        )

	emit:
		clinical_vcf  = ADD_MOST_SEVERE_CSQ_CLIN.out.vcf
		research_vcf  = ADD_MOST_SEVERE_CSQ_RES.out.vcf
        versions      = ch_versions.ifEmpty(null)     // channel: [ versions.yml ]
}
