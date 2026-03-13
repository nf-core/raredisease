//
// A subworkflow to annotate snvs
//

include { BCFTOOLS_ANNOTATE                     } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_VIEW                         } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_ANNOTATE as RENAME_CHRNAMES  } from '../../../modules/nf-core/bcftools/annotate/main'
include { CADD                                  } from '../../../modules/nf-core/cadd/main'
include { GAWK as REFERENCE_TO_CADD_CHRNAMES    } from '../../../modules/nf-core/gawk/main'
include { GAWK as CADD_TO_REFERENCE_CHRNAMES    } from '../../../modules/nf-core/gawk/main'
include { TABIX_TABIX as TABIX_CADD             } from '../../../modules/nf-core/tabix/tabix/main'

workflow ANNOTATE_CADD {

    take:
        ch_cadd_resources // channel: [mandatory] [ path(dir) ]
        ch_fai            // channel: [optional]  [ path(fai) ]
        ch_header         // channel: [mandatory] [ path(txt) ]
        ch_vcf            // channel: [mandatory] [ val(meta), path(vcfs), path(idx) ]
        val_genome        //  string: GRCh37 or GRCh38

    main:
        ch_rename_chrs    = channel.value([[]])

        if (val_genome.equals('GRCh38')) {

            REFERENCE_TO_CADD_CHRNAMES ( ch_fai , [], false )

            CADD_TO_REFERENCE_CHRNAMES ( ch_fai , [], false )

            CADD_TO_REFERENCE_CHRNAMES.out.output.map { _meta, txt -> txt }
                .set { ch_rename_chrs }

            ch_vcf
                .map { meta, vcf, tbi -> [ meta, vcf, tbi, [], [], [], [] ] }
                .combine(REFERENCE_TO_CADD_CHRNAMES.out.output.map { _meta, txt -> txt })
                .set { ch_rename_chrnames_in }

            RENAME_CHRNAMES ( ch_rename_chrnames_in )

            RENAME_CHRNAMES.out.vcf
                .map { meta, vcf -> [ meta, vcf, [] ] }
                .set { ch_vcf }
        }

        BCFTOOLS_VIEW(ch_vcf, [], [], [])

        CADD(BCFTOOLS_VIEW.out.vcf, ch_cadd_resources, [[:], []])

        TABIX_CADD(CADD.out.tsv)

        ch_vcf
            .join(CADD.out.tsv)
            .join(TABIX_CADD.out.index)
            .map { meta, vcf, tbi, ann, ann_tbi  -> [ meta, vcf, tbi, ann, ann_tbi, [] ] }
            .combine(ch_header.map { _meta, header -> header })
            .combine(ch_rename_chrs)
            .set { ch_annotate_in }

        BCFTOOLS_ANNOTATE(ch_annotate_in)

    emit:
        tbi  = BCFTOOLS_ANNOTATE.out.tbi // channel: [ val(meta), path(tbi) ]
        vcf  = BCFTOOLS_ANNOTATE.out.vcf // channel: [ val(meta), path(vcf) ]
}
