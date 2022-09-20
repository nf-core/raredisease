//
// A subworkflow to call SNVs by sentieon dnascope with a machine learning model.
//

include { ADD_MOST_SEVERE_CSQ } from '../../modules/local/add_most_severe_consequence'

workflow ANNOTATE_CSQ {
	take:
		clinical              // channel: [ val(meta), vcf ]
		research              // channel: [ val(meta), vcf ]
		variant_consequences  // path: consequences.txt

	main:
		ch_versions = Channel.empty()

        ch_input = clinical.join(research)

        ADD_MOST_SEVERE_CSQ (
            ch_input,
            variant_consequences
        )

	emit:
		vcf		 = ADD_MOST_SEVERE_CSQ.out.vcf
        versions = ch_versions.ifEmpty(null)     // channel: [ versions.yml ]
}
