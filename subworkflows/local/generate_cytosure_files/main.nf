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
        ch_sample_id_map  // channel: [optional] [val(id), val(id)]
        ch_tbi            // channel: [mandatory] [ val(meta), path(vcf_index) ]
        ch_vcf            // channel: [mandatory] [ val(meta), path(vcf) ]
        val_sample_id_map // string: path to sample_id_map file

    main:
        ch_reheader_out = channel.empty()

        TIDDIT_COV_VCF2CYTOSURE (ch_bam_bai, [[],[]])

        // Build channel: [val(sample_meta), path(vcf), path(vcf_index)]
        ch_vcf.join( ch_tbi, failOnMismatch: true )
            .set { ch_vcf_tbi }

        ch_bam_bai.combine(ch_vcf_tbi)
            .map {
                meta_sample, _bam, _bai, _meta_case, vcf, tbi ->
                def id_meta = ['id':meta_sample.sample]
                def sex_meta = ['sex':meta_sample.sex]
                return [ id_meta, sex_meta, vcf, tbi ]
            }
            .join(ch_sample_id_map, remainder: true)
            .branch { id_meta, sex_meta, vcf, tbi, samplemap  ->
                id: samplemap.equals(null)
                    return [id_meta + [custid:id_meta.id] + sex_meta, vcf, tbi]
                custid: !(samplemap.equals(null))
                    return [id_meta + [custid:samplemap] + sex_meta, vcf, tbi]
            }
            .set { ch_for_mix }

        channel.empty()
            .mix(ch_for_mix.id, ch_for_mix.custid)
            .set { ch_sample_vcf }

        // Split vcf into sample vcf:s and frequency filter
        SPLIT_AND_FILTER_SV_VCF ( ch_sample_vcf, [], [], [] )

        if (!val_sample_id_map.equals(null)) {

            SPLIT_AND_FILTER_SV_VCF.out.vcf
                .map { meta, vcf -> return [meta, vcf, [], []]}
                .set { ch_reheader_in }

            BCFTOOLS_REHEADER_SV_VCF ( ch_reheader_in, [[:],[]] ).vcf
                .set {ch_reheader_out}

        }

        SPLIT_AND_FILTER_SV_VCF.out.vcf
            .join(ch_reheader_out, remainder: true)
            .branch { meta, filteredvcf, reheaderedvcf  ->
                split: reheaderedvcf.equals(null)
                    return [meta, filteredvcf]
                reheader: !(reheaderedvcf.equals(null))
                    return [meta, reheaderedvcf]
            }
            .set { ch_for_mix }

        channel.empty()
            .mix(ch_for_mix.split, ch_for_mix.reheader)
            .toSortedList { a, b -> a[0].id <=> b[0].id }
            .flatMap()
            .set { ch_vcf2cytosure_in }

        TIDDIT_COV_VCF2CYTOSURE.out.cov
            .toSortedList { a, b -> a[0].id <=> b[0].id }
            .flatMap()
            .set { ch_cov2cytosure_in }

        VCF2CYTOSURE (
            ch_vcf2cytosure_in,
            ch_cov2cytosure_in,
            [[:], []], [[:], []],
            ch_blacklist
        )

    emit:
        cgh = VCF2CYTOSURE.out.cgh // channel: [ val(meta), path(cgh) ]
}
