//
// Convert VCF with structural variations to the “.CGH” format used by the CytoSure Interpret Software
//

include { BCFTOOLS_VIEW as SPLIT_AND_FILTER_SV_VCF      } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_SV_VCF } from '../../modules/nf-core/bcftools/reheader/main'
include { TIDDIT_COV as TIDDIT_COV_VCF2CYTOSURE         } from '../../modules/nf-core/tiddit/cov/main'
include { VCF2CYTOSURE                                  } from '../../modules/nf-core/vcf2cytosure/main'

workflow GENERATE_CYTOSURE_FILES {
    take:
        ch_vcf           // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_tbi           // channel: [mandatory] [ val(meta), path(vcf_index) ]
        ch_bam           // channel: [mandatory] [ val(meta), path(bam) ]
        ch_sample_id_map // channel: [optional] [val(id), val(id)]
        ch_blacklist     // channel: [optional] [path(blacklist)]

    main:
        ch_versions = Channel.empty()

        TIDDIT_COV_VCF2CYTOSURE (ch_bam, [[],[]])

        // Build channel: [val(sample_meta), path(vcf), path(vcf_index)]
        ch_vcf.join( ch_tbi, failOnMismatch: true )
            .set { ch_vcf_tbi }

        ch_bam.combine(ch_vcf_tbi).map {
            meta_sample, bam, meta_case, vcf, tbi ->
            new_meta = ['id':meta_sample.sample]
            return [ new_meta, vcf, tbi ]
        }.set { ch_sample_vcf }

        // Split vcf into sample vcf:s and frequency filter
        SPLIT_AND_FILTER_SV_VCF ( ch_sample_vcf, [], [], [] )

        SPLIT_AND_FILTER_SV_VCF.out.vcf
            .join(ch_sample_id_map)
            .map { meta, vcf, custid -> return [meta+[custid:custid], vcf, [], []]}
            .set { ch_reheader_in }

        BCFTOOLS_REHEADER_SV_VCF ( ch_reheader_in, [[:],[]] )

        VCF2CYTOSURE (
            BCFTOOLS_REHEADER_SV_VCF.out.vcf,
            TIDDIT_COV_VCF2CYTOSURE.out.cov,
            [[:], []], [[:], []],
            ch_blacklist
        )

        ch_versions = ch_versions.mix(TIDDIT_COV_VCF2CYTOSURE.out.versions.first())
        ch_versions = ch_versions.mix(SPLIT_AND_FILTER_SV_VCF.out.versions.first())
        ch_versions = ch_versions.mix(VCF2CYTOSURE.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_REHEADER_SV_VCF.out.versions.first())

    emit:
        versions = ch_versions // channel: [ versions.yml ]
}
