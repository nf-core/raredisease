//
// Annotate with VCFanno
//

include { VCFANNO } from '../../modules/nf-core/modules/vcfanno/main'
include { BCFTOOLS_VIEW } from '../../modules/nf-core/modules/bcftools/view/main'

workflow ANNOTATION_VCFANNO {
    take:
    main:
        ch_versions = Channel.empty()
    emit:
        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}