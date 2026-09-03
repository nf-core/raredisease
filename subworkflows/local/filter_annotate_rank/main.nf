include { VCF_FILTER_BCFTOOLS_FILTERVEP as GENERATE_CLINICAL_SET } from '../vcf_filter_bcftools_filtervep'
include { ANNOTATE_CSQ_PLI                                       } from '../annotate_consequence_pli'
include { RANK_VARIANTS                                          } from '../rank_variants'

workflow FILTER_ANNOTATE_RANK {
    take:
    ch_hgnc_ids                // channel: [ path(txt) ]
    ch_pedfile                 // channel: [ path(ped) ]
    ch_reduced_penetrance      // channel: [ path(txt) ]
    ch_score_config            // channel: [ path(ini) ]
    ch_variant_consequences    // channel: [ path(txt) ]
    ch_vcf                     // channel: [ val(meta), path(vcf) ]  - VEP-annotated vcf entering this stage
    filter_with_bcftools       // boolean
    filter_with_filter_vep     // boolean
    process_with_sort          // boolean - only meaningful when run_rank=true
    run_rank                   // boolean - false for ME (no ranking step)
    skip_generate_clinical_set // boolean
    val_type_label             // string  - used in the "skipping ranking" warn message; only meaningful when run_rank=true

    main:
    ch_clin_research_vcf = ch_vcf
        .multiMap { meta, vcf ->
            clinical: [ meta + [ set: "clinical" ], vcf ]
            research: [ meta + [ set: "research" ], vcf ]
        }

    ch_clinical_vcf = channel.empty()
    if (!skip_generate_clinical_set) {
        GENERATE_CLINICAL_SET(
            ch_clin_research_vcf.clinical,
            ch_hgnc_ids,
            filter_with_bcftools,
            filter_with_filter_vep
        )
        ch_clinical_vcf = GENERATE_CLINICAL_SET.out.vcf
    }

    ch_ann_csq_in = ch_clinical_vcf.mix(ch_clin_research_vcf.research)

    ANNOTATE_CSQ_PLI(
        ch_variant_consequences,
        ch_ann_csq_in,
        !run_rank
    )

    ch_out_vcf = ANNOTATE_CSQ_PLI.out.vcf_ann
    ch_out_tbi = ANNOTATE_CSQ_PLI.out.tbi

    if (run_rank) {
        ch_rank_in = ANNOTATE_CSQ_PLI.out.vcf_ann
            .filter { meta, _vcf ->
                if (meta.probands.size()==0) {
                    log.warn("Skipping ${val_type_label} ranking since no affected samples are detected in the case")
                }
                meta.probands.size()>0
            }

        RANK_VARIANTS(
            ch_pedfile,
            ch_reduced_penetrance,
            ch_score_config,
            ch_rank_in,
            process_with_sort
        )
        ch_out_vcf = RANK_VARIANTS.out.vcf
        ch_out_tbi = RANK_VARIANTS.out.tbi
    }

    emit:
    vcf = ch_out_vcf
    tbi = ch_out_tbi
}
