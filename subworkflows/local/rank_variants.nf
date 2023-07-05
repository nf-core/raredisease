//
// A subworkflow to score and rank variants.
//

include { GENMOD_ANNOTATE  } from '../../modules/nf-core/genmod/annotate/main'
include { GENMOD_MODELS    } from '../../modules/nf-core/genmod/models/main'
include { GENMOD_SCORE     } from '../../modules/nf-core/genmod/score/main'
include { GENMOD_COMPOUND  } from '../../modules/nf-core/genmod/compound/main'
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'

workflow RANK_VARIANTS {

    take:
        ch_vcf                // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_pedfile            // channel: [mandatory] [ path(ped) ]
        ch_reduced_penetrance // channel: [mandatory] [ path(pentrance) ]
        ch_score_config       // channel: [mandatory] [ path(ini) ]

    main:
        ch_versions = Channel.empty()

        GENMOD_ANNOTATE(ch_vcf)

        GENMOD_MODELS(GENMOD_ANNOTATE.out.vcf, ch_pedfile, ch_reduced_penetrance)

        GENMOD_SCORE(GENMOD_MODELS.out.vcf, ch_pedfile, ch_score_config)

        GENMOD_COMPOUND(GENMOD_SCORE.out.vcf)

        TABIX_BGZIPTABIX (GENMOD_COMPOUND.out.vcf)

        ch_versions = ch_versions.mix(GENMOD_ANNOTATE.out.versions)
        ch_versions = ch_versions.mix(GENMOD_MODELS.out.versions)
        ch_versions = ch_versions.mix(GENMOD_SCORE.out.versions)
        ch_versions = ch_versions.mix(GENMOD_COMPOUND.out.versions)
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    emit:
        vcf      = TABIX_BGZIPTABIX.out.gz_tbi.map { meta, vcf, tbi -> return [ meta, vcf ] }.collect() // channel: [ val(meta), path(vcf) ]
        versions = ch_versions                                                                          // channel: [ path(versions.yml) ]
}
