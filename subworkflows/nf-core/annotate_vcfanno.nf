//
// Annotate with VCFanno
//

include { VCFANNO } from '../../modules/nf-core/modules/vcfanno/main'

workflow ANNOTATE_VCFANNO {
    take:
        toml            // channel: path(toml)
        vcf             // channel: [ val(meta), path(vcf), path(tbi) ]
        resource_dir    // channel: path(resource_dir)

    main:
        ch_versions = Channel.empty()

        ch_placeholder = vcf.map { meta, vcf, idx -> vcf = []; [meta, vcf] }
        VCFANNO (vcf, ch_placeholder, toml, resource_dir)
        ch_versions = ch_versions.mix(VCFANNO.out.versions)

    emit:
        annotated_vcf          = VCFANNO.out.vcf                // channel: [ val(meta), path(*.vcf.gz) ]
        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
