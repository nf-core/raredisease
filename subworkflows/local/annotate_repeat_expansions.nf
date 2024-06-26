//
// Annotate repeat expansions
//

include { BCFTOOLS_VIEW as COMPRESS_STRANGER           } from '../../modules/nf-core/bcftools/view/main'
include { STRANGER                                     } from '../../modules/nf-core/stranger/main'
include { TABIX_TABIX as INDEX_STRANGER                } from '../../modules/nf-core/tabix/tabix/main'

workflow ANNOTATE_REPEAT_EXPANSIONS {
    take:
        ch_variant_catalog // channel: [mandatory] [ path(variant_catalog.json) ]
        ch_vcf             // channel: [mandatory] [ val(meta), path(vcf) ]

    main:
        ch_versions = Channel.empty()

        // Annotate, compress and index
        STRANGER ( ch_vcf, ch_variant_catalog )
        COMPRESS_STRANGER (
            STRANGER.out.vcf.map{ meta, vcf -> [meta, vcf, [] ]},
            [], [], []
        )
        INDEX_STRANGER ( COMPRESS_STRANGER.out.vcf )

        ch_vcf_idx = COMPRESS_STRANGER.out.vcf.join(INDEX_STRANGER.out.tbi, failOnMismatch:true, failOnDuplicate:true)

        ch_versions = ch_versions.mix(STRANGER.out.versions.first())
        ch_versions = ch_versions.mix(COMPRESS_STRANGER.out.versions.first())
        ch_versions = ch_versions.mix(INDEX_STRANGER.out.versions.first())

emit:
        vcf      = ch_vcf_idx   // channel: [ val(meta), path(vcf), path(tbi) ]
        versions = ch_versions  // channel: [ path(versions.yml) ]
}
