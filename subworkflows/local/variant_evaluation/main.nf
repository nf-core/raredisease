//
// A subworkflow to evaluate variants using truth vcfs.
//

include { BCFTOOLS_REHEADER             } from '../../../modules/nf-core/bcftools/reheader/main'
include { RTGTOOLS_VCFEVAL              } from '../../../modules/nf-core/rtgtools/vcfeval/main'
include { TABIX_TABIX as TABIX_TRUTHVCF } from '../../../modules/nf-core/tabix/tabix/main'

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

        TABIX_TRUTHVCF (BCFTOOLS_REHEADER.out.vcf)

        BCFTOOLS_REHEADER.out.vcf
            .join(TABIX_TRUTHVCF.out.index)
            .set { ch_truthvcf_tbi }

        ch_snv_vcf_tbi
            .combine(ch_truthvcf_tbi)
            .map { meta, query, qidx, meta2, truth, tidx ->
                    return [meta + [samplename: meta2.samplename] , query, qidx, truth, tidx, meta2.evaluationregions, meta2.bedregions]
            }
            .set { ch_vcfeval_in }

        RTGTOOLS_VCFEVAL ( ch_vcfeval_in, ch_sdf )

    emit:
        baseline_tbi        = RTGTOOLS_VCFEVAL.out.baseline_tbi // channel: [ val(meta), path(tbi) ]
        baseline_vcf        = RTGTOOLS_VCFEVAL.out.baseline_vcf // channel: [ val(meta), path(vcf) ]
        false_negatives_tbi = RTGTOOLS_VCFEVAL.out.fn_tbi       // channel: [ val(meta), path(tbi) ]
        false_negatives_vcf = RTGTOOLS_VCFEVAL.out.fn_vcf       // channel: [ val(meta), path(vcf) ]
        false_positives_tbi = RTGTOOLS_VCFEVAL.out.fp_tbi       // channel: [ val(meta), path(tbi) ]
        false_positives_vcf = RTGTOOLS_VCFEVAL.out.fp_vcf       // channel: [ val(meta), path(vcf) ]
        non_snp_roc         = RTGTOOLS_VCFEVAL.out.non_snp_roc  // channel: [ val(meta), path(tsv) ]
        phasing             = RTGTOOLS_VCFEVAL.out.phasing      // channel: [ val(meta), path(txt) ]
        snp_roc             = RTGTOOLS_VCFEVAL.out.snp_roc      // channel: [ val(meta), path(tsv) ]
        summary             = RTGTOOLS_VCFEVAL.out.summary      // channel: [ val(meta), path(txt) ]
        true_positives_tbi  = RTGTOOLS_VCFEVAL.out.tp_tbi       // channel: [ val(meta), path(tbi) ]
        true_positives_vcf  = RTGTOOLS_VCFEVAL.out.tp_vcf       // channel: [ val(meta), path(vcf) ]
        weighted_roc        = RTGTOOLS_VCFEVAL.out.weighted_roc // channel: [ val(meta), path(tsv) ]
}
