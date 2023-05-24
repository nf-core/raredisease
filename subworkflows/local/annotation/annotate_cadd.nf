//
// A subworkflow to annotate snvs
//

include { BCFTOOLS_ANNOTATE             } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_VIEW                 } from '../../../modules/nf-core/bcftools/view/main'
include { CADD                          } from '../../../modules/nf-core/cadd/main'
include { TABIX_TABIX as TABIX_ANNOTATE } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_CADD     } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_VIEW     } from '../../../modules/nf-core/tabix/tabix/main'

workflow ANNOTATE_CADD {

    take:
        ch_vcf            // channel: [mandatory] [ val(meta), path(vcfs) ]
        ch_index          // channel: [mandatory] [ val(meta), path(tbis) ]
        ch_header         // channel: [mandatory] [ path(txt) ]
        ch_cadd_resources // channel: [mandatory] [ path(dir) ]

    main:
        ch_versions       = Channel.empty()

        BCFTOOLS_VIEW(ch_vcf.join(ch_index), [], [], [])

        TABIX_VIEW(BCFTOOLS_VIEW.out.vcf)

        CADD(BCFTOOLS_VIEW.out.vcf, ch_cadd_resources)

        TABIX_CADD(CADD.out.tsv)

        ch_vcf
            .join(ch_index)
            .join(CADD.out.tsv)
            .join(TABIX_CADD.out.tbi)
            .combine(ch_header)
            .set { ch_annotate_in }

        BCFTOOLS_ANNOTATE(ch_annotate_in)

        TABIX_ANNOTATE (BCFTOOLS_ANNOTATE.out.vcf)

        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())
        ch_versions = ch_versions.mix(TABIX_VIEW.out.versions.first())
        ch_versions = ch_versions.mix(CADD.out.versions.first())
        ch_versions = ch_versions.mix(TABIX_CADD.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())
        ch_versions = ch_versions.mix(TABIX_ANNOTATE.out.versions.first())

    emit:
        vcf  = BCFTOOLS_ANNOTATE.out.vcf // channel: [ val(meta), path(vcf) ]
        tbi  = TABIX_ANNOTATE.out.tbi    // channel: [ val(meta), path(tbi) ]
        versions = ch_versions           // channel: [ path(versions.yml) ]
}
