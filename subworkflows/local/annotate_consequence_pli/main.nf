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
        val_publish_dir         // val: destination string, or '' to skip publishing

    main:
        CUSTOM_ADDMOSTSEVERECONSEQUENCE (ch_vcf, ch_variant_consequences)

        CUSTOM_ADDMOSTSEVEREPLI (CUSTOM_ADDMOSTSEVERECONSEQUENCE.out.vcf)

        ch_tbi_publish = channel.empty()
        if (val_index) {
            TABIX_TABIX(CUSTOM_ADDMOSTSEVEREPLI.out.vcf)
            ch_tbi_publish = TABIX_TABIX.out.index
        }

        ch_publish = channel.empty()
        if (val_publish_dir) {
            ch_publish = CUSTOM_ADDMOSTSEVEREPLI.out.vcf
                .mix(ch_tbi_publish)
                .map { meta, value -> [val_publish_dir, [meta, value]] }
        }

    emit:
        vcf_ann  = CUSTOM_ADDMOSTSEVEREPLI.out.vcf // channel: [ val(meta), path(vcf) ]
        publish  = ch_publish                       // channel: [ val(destination), val(value) ]
}
