//
// A subworkflow to add most severe consequence and pli to a vep annotated vcf
//

include { ADD_MOST_SEVERE_CSQ } from '../../modules/local/add_most_severe_consequence'
include { ADD_MOST_SEVERE_PLI } from '../../modules/local/add_most_severe_pli'

workflow ANNOTATE_CSQ_PLI {
	take:
		vcf                   // channel: [ val(meta), vcf ]
		variant_consequences  // path: consequences.txt

	main:
		ch_versions = Channel.empty()

        ADD_MOST_SEVERE_CSQ (vcf, variant_consequences)

        ADD_MOST_SEVERE_PLI (ADD_MOST_SEVERE_CSQ.out.vcf)

        ch_versions = ch_versions.mix(ADD_MOST_SEVERE_CSQ.out.versions)
        ch_versions = ch_versions.mix(ADD_MOST_SEVERE_PLI.out.versions)

	emit:
		vcf_ann  = ADD_MOST_SEVERE_PLI.out.vcf
        versions = ch_versions.ifEmpty(null)     // channel: [ versions.yml ]
}
