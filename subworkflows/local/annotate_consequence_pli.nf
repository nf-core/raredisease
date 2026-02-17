//
// A subworkflow to add most severe consequence and pli to a vep annotated vcf
//

include { CUSTOM_ADDMOSTSEVERECONSEQUENCE } from '../../modules/nf-core/custom/addmostsevereconsequence'
include { CUSTOM_ADDMOSTSEVEREPLI         } from '../../modules/nf-core/custom/addmostseverepli'
include { TABIX_TABIX                     } from '../../modules/nf-core/tabix/tabix/main'

workflow ANNOTATE_CSQ_PLI {
    take:
        ch_variant_consequences // channel: [mandatory] [ path(consequences) ]
        ch_vcf                  // channel: [mandatory] [ val(meta), path(vcf) ]
        val_index               // bool

    main:
        CUSTOM_ADDMOSTSEVERECONSEQUENCE (ch_vcf, ch_variant_consequences)

        CUSTOM_ADDMOSTSEVEREPLI (CUSTOM_ADDMOSTSEVERECONSEQUENCE.out.vcf)

        if (val_index) {
            TABIX_TABIX(CUSTOM_ADDMOSTSEVEREPLI.out.vcf)
        }

    emit:
        vcf_ann  = CUSTOM_ADDMOSTSEVEREPLI.out.vcf // channel: [ val(meta), path(vcf) ]
}
