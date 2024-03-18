//
// A subworkflow to score and rank variants.
//

include { GENMOD_ANNOTATE  } from '../../modules/nf-core/genmod/annotate/main'
include { GENMOD_MODELS    } from '../../modules/nf-core/genmod/models/main'
include { GENMOD_SCORE     } from '../../modules/nf-core/genmod/score/main'
include { GENMOD_COMPOUND  } from '../../modules/nf-core/genmod/compound/main'
include { BCFTOOLS_SORT    } from '../../modules/nf-core/bcftools/sort/main'
include { TABIX_BGZIP      } from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX      } from '../../modules/nf-core/tabix/tabix/main'

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

        BCFTOOLS_SORT(GENMOD_COMPOUND.out.vcf) // SV file needs to be sorted before indexing

        TABIX_BGZIP(GENMOD_COMPOUND.out.vcf) //run only for SNVs

        ch_vcf = TABIX_BGZIP.out.output.mix(BCFTOOLS_SORT.out.vcf)

        TABIX_TABIX (ch_vcf)

        ch_versions = ch_versions.mix(GENMOD_ANNOTATE.out.versions)
        ch_versions = ch_versions.mix(GENMOD_MODELS.out.versions)
        ch_versions = ch_versions.mix(GENMOD_SCORE.out.versions)
        ch_versions = ch_versions.mix(GENMOD_COMPOUND.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)
        ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions)
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
        vcf      = ch_vcf       // channel: [ val(meta), path(vcf) ]
        versions = ch_versions  // channel: [ path(versions.yml) ]
}
