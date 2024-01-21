//
// A subworkflow to evaluate variants using truth vcfs.
//

include { RTGTOOLS_VCFEVAL                 } from '../../modules/nf-core/rtgtools/vcfeval/main'
include { BCFTOOLS_REHEADER                } from '../../modules/nf-core/bcftools/reheader/main'
include { TABIX_TABIX as TABIX_TRUTHVCF    } from '../../modules/nf-core/tabix/tabix/main'

workflow VARIANT_EVALUATION {

    take:
        ch_snv_vcf_tbi     // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_genome_fai      // channel: [mandatory] [ val(meta), path(fai) ]
        ch_rtg_truthvcfs   // channel: [mandatory] [ val(meta), path(dbs) ]

    main:
        ch_versions = Channel.empty()

        ch_rtg_truthvcfs
            .splitCsv ( header:true )
            .map { row ->
                return [[id:row.samplename[0], bed:row.bed[0]], row.vcf[0], []]
            }
            .set { ch_rtgvcfs_dbs }

        BCFTOOLS_REHEADER (ch_rtgvcfs_dbs, [[:],[]])

        TABIX_TRUTHVCF (BCFTOOLS_REHEADER.out.vcf)

        BCFTOOLS_REHEADER.out.vcf
            .join(TABIX_TRUTHVCF.out.tbi)
            .map { meta, vcf, tbi -> return [vcf, tbi, meta.bed]}
            .set { ch_truthvcf_tbi }

        ch_snv_vcf_tbi
            .combine(ch_truthvcf_tbi)
            .map { meta, query, qidx, truth, tidx, tbed ->
                    return [meta, query, qidx, truth, tidx, tbed, []]
            }
            .set { ch_vcfeval_in }

        RTGTOOLS_VCFEVAL ( ch_vcfeval_in )

        // ch_versions = ch_versions.mix(BUILD_BED.out.versions)
        // ch_versions = ch_versions.mix(GATK4_SPLITINTERVALS.out.versions)

    emit:
        versions        = ch_versions                   // channel: [ path(versions.yml) ]
}
