//
// A subworkflow to add pli to a vep annotated vcf
//

include { CUSTOM_ADDMOSTSEVEREPLI } from '../../../modules/local/custom/addmostseverepli'
include { TABIX_TABIX             } from '../../../modules/nf-core/tabix/tabix/main'

workflow ANNOTATE_PLI {
    take:
        ch_vcf    // channel: [mandatory] [ val(meta), path(vcf) ]
        val_index // bool

    main:
        CUSTOM_ADDMOSTSEVEREPLI (ch_vcf)

        ch_tbi = channel.empty()
        if (val_index) {
            TABIX_TABIX(CUSTOM_ADDMOSTSEVEREPLI.out.vcf)
            ch_tbi = TABIX_TABIX.out.index
        }

    emit:
        tbi     = ch_tbi                           // channel: [ val(meta), path(tbi) ]
        vcf_ann = CUSTOM_ADDMOSTSEVEREPLI.out.vcf  // channel: [ val(meta), path(vcf) ]
}
