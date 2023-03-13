//
// A subworkflow to score and rank variants.
//

include { GENMOD_ANNOTATE  } from '../../modules/nf-core/genmod/annotate/main'
include { GENMOD_MODELS    } from '../../modules/nf-core/genmod/models/main'
include { GENMOD_SCORE     } from '../../modules/nf-core/genmod/score/main'
include { GENMOD_COMPOUND  } from '../../modules/nf-core/genmod/compound/main'

workflow RANK_VARIANTS {

    take:
        ch_vcf                // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_ped                // channel: [mandatory] [ path(ped) ]
        ch_reduced_penetrance // channel: [mandatory] [ path(pentrance) ]
        ch_score_config       // channel: [mandatory] [ path(ini) ]

    main:
        ch_versions = Channel.empty()

        GENMOD_ANNOTATE(ch_vcf)

        GENMOD_MODELS(GENMOD_ANNOTATE.out.vcf, ch_ped, ch_reduced_penetrance)

        GENMOD_SCORE(GENMOD_MODELS.out.vcf, ch_ped, ch_score_config)

        GENMOD_COMPOUND(GENMOD_SCORE.out.vcf)

        ch_versions = ch_versions.mix(GENMOD_ANNOTATE.out.versions)
        ch_versions = ch_versions.mix(GENMOD_MODELS.out.versions)
        ch_versions = ch_versions.mix(GENMOD_SCORE.out.versions)
        ch_versions = ch_versions.mix(GENMOD_COMPOUND.out.versions)

    emit:
        vcf         = GENMOD_COMPOUND.out.vcf  // channel: [ val(meta), path(vcf) ]
        versions    = ch_versions              // channel: [ path(versions.yml) ]
}
