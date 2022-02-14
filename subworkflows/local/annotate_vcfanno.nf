//
// Annotate with VCFanno
//

include { VCFANNO } from '../../modules/local/vcfanno/main'
include { BCFTOOLS_ANNOTATE } from '../../modules/local/bcftools_annotate/main'

workflow ANNOTATE_VCFANNO {
    take:
        vcf             // channel: [ val(meta), path(vcf), path(tbi) ]
        resource_dir    // channel: path(resource_dir)

    main:
        ch_versions = Channel.empty()

        VCFANNO (vcf, resource_dir)
        BCFTOOLS_ANNOTATE ( VCFANNO.out.vcf )
        ch_versions = ch_versions.mix(VCFANNO.out.versions)

    emit:
        annotated_vcf          = BCFTOOLS_ANNOTATE.out.vcf      // channel: [ val(meta), path(*.vcf.gz) ]
        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}