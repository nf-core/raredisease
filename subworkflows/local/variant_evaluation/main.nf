//
// A subworkflow to evaluate variants using truth vcfs.
//

include { BCFTOOLS_REHEADER             } from '../../../modules/nf-core/bcftools/reheader/main'
include { RTGTOOLS_VCFEVAL              } from '../../../modules/nf-core/rtgtools/vcfeval/main'

workflow VARIANT_EVALUATION {

    take:
        ch_rtg_truthvcfs   // channel: [mandatory] [ val(meta), path(dbs) ]
        ch_sdf             // channel: [mandatory] [ val(meta), path(sdf) ]
        ch_snv_vcf_tbi     // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ]

    main:
        ch_rtg_truthvcfs
            .splitCsv ( header:true )
            .map { row ->
                def evregions  = row.evaluationregions[0].isEmpty() ? [] : row.evaluationregions[0]
                def bedregions = row.bedregions[0].isEmpty()        ? [] : row.bedregions[0]
                return [[samplename:row.samplename[0], bedregions:bedregions, evaluationregions:evregions], row.vcf[0], [], []]
            }
            .set { ch_rtgvcfs_dbs }

        BCFTOOLS_REHEADER (ch_rtgvcfs_dbs, [[:],[]])

        BCFTOOLS_REHEADER.out.vcf
            .join(BCFTOOLS_REHEADER.out.index, failOnMismatch:true, failOnDuplicate:true)
            .set { ch_truthvcf_tbi }

        ch_snv_vcf_tbi
            .combine(ch_truthvcf_tbi)
            .map { meta, query, qidx, meta2, truth, tidx ->
                    return [meta + [samplename: meta2.samplename] , query, qidx, truth, tidx, meta2.evaluationregions, meta2.bedregions]
            }
            .set { ch_vcfeval_in }

        RTGTOOLS_VCFEVAL ( ch_vcfeval_in, ch_sdf )

        ch_publish = RTGTOOLS_VCFEVAL.out.tp_vcf
            .mix(RTGTOOLS_VCFEVAL.out.tp_tbi)
            .mix(RTGTOOLS_VCFEVAL.out.fn_vcf)
            .mix(RTGTOOLS_VCFEVAL.out.fn_tbi)
            .mix(RTGTOOLS_VCFEVAL.out.fp_vcf)
            .mix(RTGTOOLS_VCFEVAL.out.fp_tbi)
            .mix(RTGTOOLS_VCFEVAL.out.baseline_vcf)
            .mix(RTGTOOLS_VCFEVAL.out.baseline_tbi)
            .mix(RTGTOOLS_VCFEVAL.out.snp_roc)
            .mix(RTGTOOLS_VCFEVAL.out.non_snp_roc)
            .mix(RTGTOOLS_VCFEVAL.out.weighted_roc)
            .mix(RTGTOOLS_VCFEVAL.out.summary)
            .mix(RTGTOOLS_VCFEVAL.out.phasing)
            .map { meta, value -> ['rtgvcfeval/', [meta, value]] }

    emit:
        publish = ch_publish // channel: [ val(destination), val(value) ]
}
