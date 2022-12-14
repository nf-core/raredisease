//
// A subworkflow to score and rank variants.
//

include { GENMOD_ANNOTATE  } from '../../modules/nf-core/genmod/annotate/main'
include { GENMOD_MODELS    } from '../../modules/nf-core/genmod/models/main'
include { GENMOD_SCORE     } from '../../modules/nf-core/genmod/score/main'
include { GENMOD_COMPOUND  } from '../../modules/nf-core/genmod/compound/main'

workflow RANK_VARIANTS {

    take:
        vcf                // channel: [ val(meta), path(vcf) ]
        ped
        reduced_penetrance
        score_config

    main:
        ch_versions = Channel.empty()

        GENMOD_ANNOTATE(vcf)
        ch_versions = ch_versions.mix(GENMOD_ANNOTATE.out.versions)
        GENMOD_MODELS(GENMOD_ANNOTATE.out.vcf, ped, reduced_penetrance)
        ch_versions = ch_versions.mix(GENMOD_MODELS.out.versions)
        GENMOD_SCORE(GENMOD_MODELS.out.vcf, ped, score_config)
        ch_versions = ch_versions.mix(GENMOD_SCORE.out.versions)
        GENMOD_COMPOUND(GENMOD_SCORE.out.vcf)
        ch_versions = ch_versions.mix(GENMOD_COMPOUND.out.versions)

    emit:
        vcf                    = GENMOD_COMPOUND.out.vcf
        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
