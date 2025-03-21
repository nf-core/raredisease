//
// A subworkflow to evaluate variants using truth vcfs.
//

include { RTGTOOLS_VCFEVAL                 } from '../../modules/nf-core/rtgtools/vcfeval/main'
include { BCFTOOLS_REHEADER                } from '../../modules/nf-core/bcftools/reheader/main'
include { TABIX_TABIX as TABIX_TRUTHVCF    } from '../../modules/nf-core/tabix/tabix/main'

workflow VARIANT_EVALUATION {

    take:
        ch_snv_vcf_tbi     // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ]
        ch_genome_fai      // channel: [mandatory] [ val(meta), path(fai) ]
        ch_rtg_truthvcfs   // channel: [mandatory] [ val(meta), path(dbs) ]
        ch_sdf             // channel: [mandatory] [ val(meta), path(sdf) ]

    main:
        ch_versions = Channel.empty()

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
            .join(TABIX_TRUTHVCF.out.tbi)
            .set { ch_truthvcf_tbi }

        ch_snv_vcf_tbi
            .combine(ch_truthvcf_tbi)
            .map { meta, query, qidx, meta2, truth, tidx ->
                    return [meta + [samplename: meta2.samplename] , query, qidx, truth, tidx, meta2.evaluationregions, meta2.bedregions]
            }
            .set { ch_vcfeval_in }

        RTGTOOLS_VCFEVAL ( ch_vcfeval_in, ch_sdf )

        ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions)
        ch_versions = ch_versions.mix(TABIX_TRUTHVCF.out.versions)
        ch_versions = ch_versions.mix(RTGTOOLS_VCFEVAL.out.versions)

    emit:
        versions        = ch_versions // channel: [ path(versions.yml) ]
}
