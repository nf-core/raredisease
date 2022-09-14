//
// A subworkflow to call SNVs by sentieon dnascope with a machine learning model.
//

include { VCFPARSER } from '../../modules/local/vcfparser'

workflow VCFPARSER_CSQ {
	take:
		clinical              // channel: [ val(meta), vcf ]
		research              // channel: [ val(meta), vcf ]
		variant_consequences  // path: consequences.txt

	main:
		ch_versions = Channel.empty()

        ch_input = clinical.join(research)

        VCFPARSER (
            ch_input,
            variant_consequences
        )

	emit:
		vcf		 = VCFPARSER.out.vcf
        versions = ch_versions.ifEmpty(null)     // channel: [ versions.yml ]
}
