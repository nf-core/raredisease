//
// A subworkflow to add most severe consequence and pli to a vep annotated vcf
//

include { CUSTOM_ADDMOSTSEVERECONSEQUENCE } from '../../../modules/local/custom/addmostsevereconsequence'
include { CUSTOM_ADDMOSTSEVEREPLI         } from '../../../modules/local/custom/addmostseverepli'
include { TABIX_TABIX                     } from '../../../modules/nf-core/tabix/tabix/main'

workflow ANNOTATE_CSQ_PLI {
    take:
        ch_variant_consequences // channel: [mandatory] [ path(consequences) ]
        ch_vcf                  // channel: [mandatory] [ val(meta), path(vcf) ]
        val_index               // bool

    main:
        CUSTOM_ADDMOSTSEVERECONSEQUENCE (ch_vcf, ch_variant_consequences)

        CUSTOM_ADDMOSTSEVEREPLI (CUSTOM_ADDMOSTSEVERECONSEQUENCE.out.vcf)

        ch_tbi = channel.empty()
        if (val_index) {
            TABIX_TABIX(CUSTOM_ADDMOSTSEVEREPLI.out.vcf)
            ch_tbi = TABIX_TABIX.out.index
        }

    emit:
        tbi     = ch_tbi                           // channel: [ val(meta), path(tbi) ]
        vcf_ann = CUSTOM_ADDMOSTSEVEREPLI.out.vcf  // channel: [ val(meta), path(vcf) ]
}
