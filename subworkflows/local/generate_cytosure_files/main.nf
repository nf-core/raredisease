//
// Convert VCF with structural variations to the “.CGH” format used by the CytoSure Interpret Software
//

include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_SV_VCF } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_VIEW as SPLIT_AND_FILTER_SV_VCF      } from '../../../modules/nf-core/bcftools/view/main'
include { TIDDIT_COV as TIDDIT_COV_VCF2CYTOSURE         } from '../../../modules/nf-core/tiddit/cov/main'
include { VCF2CYTOSURE                                  } from '../../../modules/nf-core/vcf2cytosure/main'

workflow GENERATE_CYTOSURE_FILES {
    take:
        ch_bam_bai        // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_blacklist      // channel: [optional] [path(blacklist)]
        ch_tbi            // channel: [mandatory] [ val(meta), path(vcf_index) ]
        ch_vcf            // channel: [mandatory] [ val(meta), path(vcf) ]

    main:
        TIDDIT_COV_VCF2CYTOSURE (ch_bam_bai, [[],[]])

        // Build channel: [val(sample_meta), path(vcf), path(vcf_index)]
        ch_vcf_tbi = ch_vcf.join( ch_tbi, failOnMismatch: true )

        ch_sample_vcf = ch_bam_bai.combine(ch_vcf_tbi)
            .map {
                meta_sample, _bam, _bai, _meta_case, vcf, tbi ->
                def new_meta = ['id':meta_sample.sample, 'sex':meta_sample.sex, 'custid':meta_sample.customer_id ?: meta_sample.sample]
                return [ new_meta, vcf, tbi ]
            }

        // Split vcf into sample vcf:s and frequency filter
        SPLIT_AND_FILTER_SV_VCF ( ch_sample_vcf, [], [], [] )

        // Only rows with a distinct customer id need their VCF header renamed
        ch_split = SPLIT_AND_FILTER_SV_VCF.out.vcf
            .branch { meta, vcf ->
                reheader: meta.custid != meta.id
                as_is:    true
            }

        ch_reheader_in = ch_split.reheader
            .map { meta, vcf -> return [meta, vcf, [], []]}

        BCFTOOLS_REHEADER_SV_VCF ( ch_reheader_in, [[:],[]] )

        ch_vcf2cytosure_in = ch_split.as_is
            .mix(BCFTOOLS_REHEADER_SV_VCF.out.vcf)
            .toSortedList { a, b -> a[0].id <=> b[0].id }
            .flatMap()

        ch_cov2cytosure_in = TIDDIT_COV_VCF2CYTOSURE.out.cov
            .toSortedList { a, b -> a[0].id <=> b[0].id }
            .flatMap()

        VCF2CYTOSURE (
            ch_vcf2cytosure_in,
            ch_cov2cytosure_in,
            [[:], []], [[:], []],
            ch_blacklist
        )

    emit:
        cgh = VCF2CYTOSURE.out.cgh // channel: [ val(meta), path(cgh) ]
}
