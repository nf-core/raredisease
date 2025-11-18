//
// A subworkflow to score and rank variants.
//

include { GENMOD_ANNOTATE                                       } from '../../../modules/nf-core/genmod/annotate/main'
include { GENMOD_MODELS                                         } from '../../../modules/nf-core/genmod/models/main'
include { GENMOD_SCORE                                          } from '../../../modules/nf-core/genmod/score/main'
include { GENMOD_SCORE as GENMOD_SCORE_FOR_GICAM                } from '../../../modules/nf-core/genmod/score/main'
include { GENMOD_COMPOUND                                       } from '../../../modules/nf-core/genmod/compound/main'
include { MIVMIR_INFER                                          } from '../../../modules/local/mivmir/main'
include { GICAM_INFER                                           } from '../../../modules/local/gicam/main'
include { BCFTOOLS_SORT                                         } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_BGZIPTABIX                                      } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_GICAM            } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_GENMOD_GICAM     } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_ANNOTATE as BCFTOOLS_MERGE_GENMOD_GICAM      } from '../../../modules/nf-core/bcftools/annotate/main'

workflow RANK_VARIANTS {

    take:
        ch_pedfile                    // channel: [mandatory] [ path(ped) ]
        ch_reduced_penetrance         // channel: [mandatory] [ path(pentrance) ]
        ch_score_config               // channel: [mandatory] [ path(ini) ]
        ch_vcf                        // channel: [mandatory] [ val(meta), path(vcf) ]
        process_with_sort             // Boolean
        rank_with_mivmir_gicam        // Boolean
        ch_genmod_gicam_score_config  // channel: [mandatory if rank_with_mivmir_gicam] [ path(ini) ]

    main:
        GENMOD_ANNOTATE(ch_vcf)

        ch_models_in = GENMOD_ANNOTATE.out.vcf.combine(ch_pedfile)

        GENMOD_MODELS(ch_models_in, ch_reduced_penetrance)

        ch_score_in = GENMOD_MODELS.out.vcf.combine(ch_pedfile)

        GENMOD_SCORE(ch_score_in, ch_score_config)

        // Run MIVMIR - GICAM scoring (not supported for MT SNVs and SVs)
        if (rank_with_mivmir_gicam) {
            GENMOD_SCORE_FOR_GICAM(ch_score_in, ch_genmod_gicam_score_config)
            MIVMIR_INFER(GENMOD_SCORE_FOR_GICAM.out.vcf)
            GICAM_INFER(MIVMIR_INFER.out.vcf)
            TABIX_BGZIPTABIX_GICAM(GICAM_INFER.out.vcf)
        }

        GENMOD_COMPOUND(GENMOD_SCORE.out.vcf)

        ch_sort_publish  = channel.empty()
        ch_tabix_publish = channel.empty()

        if (process_with_sort) {
            ch_vcf = BCFTOOLS_SORT(GENMOD_COMPOUND.out.vcf).vcf // SV file needs to be sorted before indexing
            ch_sort_publish = BCFTOOLS_SORT.out.vcf
                .mix(BCFTOOLS_SORT.out.tbi)
                .map { meta, value -> ['rank_and_filter/', [meta, value]] }
        } else {
            ch_vcf = TABIX_BGZIPTABIX(GENMOD_COMPOUND.out.vcf).gz_index.map {meta, vcf, _tbi -> return [meta, vcf]} //run only for SNVs
            ch_tabix_publish = TABIX_BGZIPTABIX.out.gz_index
                .map { meta, gz, tbi -> ['rank_and_filter/', [meta, gz, tbi]] }
        }

        ch_publish = ch_sort_publish.mix(ch_tabix_publish)

        // Merge Genmod and MIVMIR-GICAM scores
        if (rank_with_mivmir_gicam) {
            ch_vcf.join(TABIX_BGZIPTABIX_GICAM.out.gz_index, failOnMismatch: true)
            .map {meta, vcf_genmod, vcf_gicam, vcf_index_gicam -> return [ meta, vcf_genmod, [], vcf_gicam, vcf_index_gicam, [], [], [] ]}
            .set {ch_merge_genmod_gicam}
            BCFTOOLS_MERGE_GENMOD_GICAM(ch_merge_genmod_gicam)
            TABIX_BGZIPTABIX_GENMOD_GICAM(BCFTOOLS_MERGE_GENMOD_GICAM.out.vcf)
            TABIX_BGZIPTABIX_GENMOD_GICAM.out.gz_index.map {meta, vcf, _tbi -> return [meta, vcf]}.set {ch_vcf}
            TABIX_BGZIPTABIX_GENMOD_GICAM.out.gz_index
                .map { meta, gz, tbi -> ['rank_and_filter/', [meta, gz, tbi]] }.set {ch_publish}
        }

    emit:
        publish  = ch_publish   // channel: [ val(destination), val(value) ]
        vcf      = ch_vcf       // channel: [ val(meta), path(vcf) ]
}
