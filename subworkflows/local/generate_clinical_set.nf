//
// Generarte clinical set of variants
//

include { ENSEMBLVEP_FILTERVEP } from '../../modules/nf-core/ensemblvep/filtervep'
include { TABIX_BGZIP          } from '../../modules/nf-core/tabix/bgzip'
include { TABIX_TABIX          } from '../../modules/nf-core/tabix/tabix'

workflow GENERATE_CLINICAL_SET {
    take:
        ch_vcf      // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_hgnc_ids // channel: [mandatory] [ val(hgnc_ids) ]

    main:
        ch_versions = Channel.empty()

        ch_vcf
            .multiMap { meta, vcf ->
                clinical: [ meta + [ set: "clinical" ], vcf ]
                research: [ meta + [ set: "research" ], vcf ]
            }
            .set { ch_clin_research_vcf }

        ENSEMBLVEP_FILTERVEP(
            ch_clin_research_vcf.clinical,
            ch_hgnc_ids
        )
        .output
        .set { ch_filtervep_out }

        TABIX_BGZIP( ch_filtervep_out )

        ch_clin_research_vcf.research
            .mix( TABIX_BGZIP.out.output )
            .set { ch_clin_research_split }

        TABIX_TABIX( ch_clin_research_split )

        ch_versions = ch_versions.mix( ENSEMBLVEP_FILTERVEP.out.versions )
        ch_versions = ch_versions.mix( TABIX_BGZIP.out.versions )
        ch_versions = ch_versions.mix( TABIX_TABIX.out.versions )

    emit:
        vcf      = ch_clin_research_split // channel: [ val(meta), path(vcf) ]
        tbi      = TABIX_TABIX.out.tbi    // channel: [ val(meta), path(tbi) ]
        versions = ch_versions            // channel: [ path(versions.yml) ]
}
