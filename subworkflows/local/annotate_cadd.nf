//
// A subworkflow to annotate snvs
//

include { BCFTOOLS_VIEW } from '../../modules/nf-core/bcftools/view/main'
include { CADD          } from '../../modules/local/cadd'

workflow ANNOTATE_CADD {

    take:
        ch_vcf         // channel: [mandatory] [ val(meta), path(vcfs) ]
        ch_index       // channel: [mandatory] [ val(meta), path(tbis) ]
        ch_cadd_scores // channel: [mandatory] [ val(meta), path(dir) ]

    main:
        ch_versions       = Channel.empty()

        BCFTOOLS_VIEW(ch_vcf.join(ch_index),[],[],[])

        CADD(BCFTOOLS_VIEW.out.vcf, ch_cadd_scores)

    emit:
        vcf_ann  = Channel.empty()   // channel: [ val(meta), path(vcf) ]
        // tbi      = ch_vep_index // channel: [ val(meta), path(tbi) ]
        // versions = ch_versions  // channel: [ path(versions.yml) ]
}
