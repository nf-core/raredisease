//
// Annotate with VCFanno
//

include { VCFANNO } from '../../modules/nf-core/modules/vcfanno/main'
include { BCFTOOLS_VIEW } from '../../modules/nf-core/modules/bcftools/view/main'

workflow ANNOTATION_VCFANNO {
    take:
        vcf             // channel: [ val(meta), path(vcf), path(tbi) ]
        toml_config     // channel: path(.toml)

    main:
        ch_versions = Channel.empty()

    emit:
        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}