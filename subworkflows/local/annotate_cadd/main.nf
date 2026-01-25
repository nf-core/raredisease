//
// A subworkflow to annotate snvs
//

include { BCFTOOLS_ANNOTATE                     } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_VIEW                         } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_ANNOTATE as RENAME_CHRNAMES  } from '../../../modules/nf-core/bcftools/annotate/main'
include { CADD                                  } from '../../../modules/nf-core/cadd/main'
include { GAWK as REFERENCE_TO_CADD_CHRNAMES    } from '../../../modules/nf-core/gawk/main'
include { GAWK as CADD_TO_REFERENCE_CHRNAMES    } from '../../../modules/nf-core/gawk/main'
include { TABIX_TABIX as TABIX_ANNOTATE         } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_CADD             } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_VIEW             } from '../../../modules/nf-core/tabix/tabix/main'

workflow ANNOTATE_CADD {

    take:
        ch_cadd_resources // channel: [mandatory] [ path(dir) ]
        ch_fai            // channel: [optional]  [ path(fai) ]
        ch_header         // channel: [mandatory] [ path(txt) ]
        ch_vcf            // channel: [mandatory] [ val(meta), path(vcfs), path(idx) ]
        val_genome        //  string: GRCh37 or GRCh37

    main:
        ch_versions       = channel.empty()
        ch_rename_chrs    = channel.empty()

        if (val_genome.equals('GRCh38')) {

            REFERENCE_TO_CADD_CHRNAMES ( ch_fai , [], false )
            ch_versions = ch_versions.mix(REFERENCE_TO_CADD_CHRNAMES.out.versions)

            CADD_TO_REFERENCE_CHRNAMES ( ch_fai , [], false )
            ch_versions = ch_versions.mix(CADD_TO_REFERENCE_CHRNAMES.out.versions)

            CADD_TO_REFERENCE_CHRNAMES.out.output.map { _meta, txt -> txt }
                .set { ch_rename_chrs }

            ch_vcf
                .map { meta, vcf, tbi -> [ meta, vcf, tbi, [], [] ] }
                .set { rename_chrnames_in }

            RENAME_CHRNAMES (
                rename_chrnames_in,
                [],
                [],
                REFERENCE_TO_CADD_CHRNAMES.out.output.map { _meta, txt -> txt }
            )

            RENAME_CHRNAMES.out.vcf
                .map { meta, vcf -> [ meta, vcf, [] ] }
                .set { ch_vcf }
        }


        BCFTOOLS_VIEW(ch_vcf, [], [], [])

        TABIX_VIEW(BCFTOOLS_VIEW.out.vcf)

        CADD(BCFTOOLS_VIEW.out.vcf, ch_cadd_resources)

        TABIX_CADD(CADD.out.tsv)

        ch_vcf
            .join(CADD.out.tsv)
            .join(TABIX_CADD.out.tbi)
            .map { meta, vcf, annotations, annotations_index -> [ meta, vcf, [], annotations, annotations_index ] }
            .set { ch_annotate_in }

        BCFTOOLS_ANNOTATE(ch_annotate_in, [], ch_header.map { _meta, header -> header }, ch_rename_chrs)

        TABIX_ANNOTATE (BCFTOOLS_ANNOTATE.out.vcf)

        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)
        ch_versions = ch_versions.mix(TABIX_VIEW.out.versions)
        ch_versions = ch_versions.mix(CADD.out.versions)
        ch_versions = ch_versions.mix(TABIX_CADD.out.versions)
        ch_versions = ch_versions.mix(TABIX_ANNOTATE.out.versions)

    emit:
        tbi  = TABIX_ANNOTATE.out.tbi    // channel: [ val(meta), path(tbi) ]
        vcf  = BCFTOOLS_ANNOTATE.out.vcf // channel: [ val(meta), path(vcf) ]
        versions = ch_versions           // channel: [ path(versions.yml) ]
}
