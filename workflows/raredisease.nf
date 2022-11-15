/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRaredisease.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.bwa_index,
    params.bwamem2_index,
    params.call_interval,
    params.fasta,
    params.fasta_fai,
    params.gens_gnomad_pos,
    params.gens_interval_list,
    params.gens_pon,
    params.gnomad_af,
    params.gnomad_af_idx,
    params.gnomad_vcf,
    params.input,
    params.known_dbsnp,
    params.known_dbsnp_tbi,
    params.known_indels,
    params.known_mills,
    params.ml_model,
    params.mt_backchain_shift,
    params.mt_bwamem2_index_shift,
    params.mt_fasta_shift,
    params.mt_fai_shift,
    params.mt_intervals,
    params.mt_intervals_shift,
    params.mt_sequence_dictionary_shift,
    params.multiqc_config,
    params.reduced_penetrance,
    params.score_config_snv,
    params.score_config_sv,
    params.sequence_dictionary,
    params.target_bed,
    params.svdb_query_dbs,
    params.variant_catalog,
    params.vep_filters,
    params.vcfanno_resources,
    params.vcfanno_toml,
    params.vep_cache
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config              = params.multiqc_config              ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo                       = params.multiqc_logo                ? Channel.fromPath( params.multiqc_logo, checkIfExists: true )   : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true)  : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: local modules
//

include { FILTER_VEP as FILTER_VEP_SNV          } from '../modules/local/filter_vep'
include { FILTER_VEP as FILTER_VEP_SV           } from '../modules/local/filter_vep'
include { MAKE_PED                              } from '../modules/local/create_pedfile'

//
// MODULE: Installed directly from nf-core/modules
//

include { CUSTOM_DUMPSOFTWAREVERSIONS           } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTQC                                } from '../modules/nf-core/fastqc/main'
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'

//
// SUBWORKFLOWS
//

include { ALIGN                                 } from '../subworkflows/local/align'
include { ANALYSE_MT                            } from '../subworkflows/local/analyse_MT'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SNV   } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SV    } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_SNVS                         } from '../subworkflows/local/annotate_snvs'
include { ANNOTATE_STRUCTURAL_VARIANTS          } from '../subworkflows/local/annotate_structural_variants'
include { CALL_REPEAT_EXPANSIONS                } from '../subworkflows/local/call_repeat_expansions'
include { CALL_SNV                              } from '../subworkflows/local/call_snv'
include { CALL_STRUCTURAL_VARIANTS              } from '../subworkflows/local/call_structural_variants'
include { CHECK_INPUT                           } from '../subworkflows/local/check_input'
include { GENS                                  } from '../subworkflows/local/gens'
include { PREPARE_REFERENCES                    } from '../subworkflows/local/prepare_references'
include { QC_BAM                                } from '../subworkflows/local/qc_bam'
include { RANK_VARIANTS as RANK_VARIANTS_SNV    } from '../subworkflows/local/rank_variants'
include { RANK_VARIANTS as RANK_VARIANTS_SV     } from '../subworkflows/local/rank_variants'
include { SCATTER_GENOME                        } from '../subworkflows/local/scatter_genome'
//
// SUBWORKFLOW: Consists entirely of nf-core/modules
//


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RAREDISEASE {

    ch_versions = Channel.empty()

    // Initialize input channels
    if (params.input) {
        ch_input = Channel.fromPath(params.input)
        CHECK_INPUT (ch_input)
    } else {
        exit 1, 'Input samplesheet not specified!'
    }

    // Initialize all file channels including unprocessed vcf, bed and tab files
    ch_call_interval                  = params.call_interval          ? Channel.fromPath(params.call_interval).collect()
                                                                      : Channel.value([])
    ch_genome_fasta_no_meta           = params.fasta                  ? Channel.fromPath(params.fasta).collect()
                                                                      : ( exit 1, 'Genome fasta not specified!' )
    ch_genome_fasta_meta              = ch_genome_fasta_no_meta.map { it -> [[id:it[0].simpleName], it] }
    ch_gnomad_af_tab                  = params.gnomad_af              ? Channel.fromPath(params.gnomad_af).map{ it -> [[id:it[0].simpleName], it] }.collect()
                                                                      : Channel.value([[],[]])
    ch_gnomad_vcf_unprocessed         = params.gnomad_vcf             ? Channel.fromPath(params.gnomad_vcf).collect()
                                                                      : Channel.empty()
    ch_known_dbsnp                    = params.known_dbsnp            ? Channel.fromPath(params.known_dbsnp).map{ it -> [[id:it[0].simpleName], it] }.collect()
                                                                      : Channel.value([[],[]])
    ch_ml_model                       = params.ml_model               ? Channel.fromPath(params.ml_model).collect()
                                                                      : Channel.value([])
    ch_mt_backchain_shift             = params.mt_backchain_shift     ? Channel.fromPath(params.mt_backchain_shift).collect()
                                                                      : Channel.value([])
    ch_mt_fasta_shift_no_meta         = params.mt_fasta_shift         ? Channel.fromPath(params.mt_fasta_shift).collect()
                                                                      : Channel.empty()
    ch_mt_fasta_shift_meta            = params.mt_fasta_shift         ? ch_mt_fasta_shift_no_meta.map { it -> [[id:it[0].simpleName], it] }.collect()
                                                                      : Channel.empty()
    ch_mt_intervals                   = params.mt_intervals           ? Channel.fromPath(params.mt_intervals).collect()
                                                                      : Channel.value([])
    ch_mt_intervals_shift             = params.mt_intervals_shift     ? Channel.fromPath(params.mt_intervals_shift).collect()
                                                                      : Channel.value([])
    ch_reduced_penetrance             = params.reduced_penetrance     ? Channel.fromPath(params.reduced_penetrance).collect()
                                                                      : Channel.value([])
    ch_score_config_snv               = params.score_config_snv       ? Channel.fromPath(params.score_config_snv).collect()
                                                                      : Channel.value([])
    ch_score_config_sv                = params.score_config_sv        ? Channel.fromPath(params.score_config_sv).collect()
                                                                      : Channel.value([])
    ch_target_bed_unprocessed         = params.target_bed             ? Channel.fromPath(params.target_bed).map{ it -> [[id:it[0].simpleName], it] }.collect()
                                                                      : Channel.value([])
    ch_variant_catalog                = params.variant_catalog        ? Channel.fromPath(params.variant_catalog).collect()
                                                                      : Channel.value([])
    ch_variant_consequences           = Channel.fromPath("$projectDir/assets/variant_consequences_v1.txt", checkIfExists: true).collect()

    ch_vcfanno_resources_unprocessed  = params.vcfanno_resources      ? Channel.fromPath(params.vcfanno_resources).map{ it -> [[id:it[0].simpleName], it] }.collect()
                                                                      : Channel.value([])
    ch_vep_cache                      = params.vep_cache              ? Channel.fromPath(params.vep_cache).collect()
                                                                      : Channel.value([])
    ch_vep_filters                    = params.vep_filters            ? Channel.fromPath(params.vep_filters).collect()
                                                                      : Channel.value([])

    // Generate pedigree file
    MAKE_PED (CHECK_INPUT.out.samples.toList())
    ch_versions = ch_versions.mix(MAKE_PED.out.versions)

    // Input QC
    FASTQC (CHECK_INPUT.out.reads)
    ch_versions = ch_versions.mix(CHECK_INPUT.out.versions)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // Prepare references and indices.
    PREPARE_REFERENCES (
        ch_genome_fasta_no_meta,
        ch_genome_fasta_meta,
        params.fasta_fai,
        ch_mt_fasta_shift_no_meta,
        ch_mt_fasta_shift_meta,
        ch_gnomad_af_tab,
        ch_gnomad_vcf_unprocessed,
        ch_known_dbsnp,
        ch_target_bed_unprocessed,
        ch_vcfanno_resources_unprocessed
    )
    .set { ch_references }

    // Gather built indices or get them from the params
    ch_bait_intervals               = ch_references.bait_intervals
    ch_bwa_index                    = params.bwa_index                     ? Channel.fromPath(params.bwa_index).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                           : ( ch_references.bwa_index                ?: Channel.empty() )
    ch_bwamem2_index                = params.bwamem2_index                 ? Channel.fromPath(params.bwamem2_index).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                           : ( ch_references.bwamem2_index            ?: Channel.empty() )
    ch_bwamem2_index_mt_shift       = params.mt_bwamem2_index_shift        ? Channel.fromPath(params.mt_bwamem2_index_shift).collect()
                                                                           : ( ch_references.bwamem2_index_mt_shift   ?: Channel.empty() )
    ch_chrom_sizes                  = ch_references.chrom_sizes
    ch_genome_fai_no_meta           = params.fasta_fai                     ? Channel.fromPath(params.fasta_fai).collect()
                                                                           : ( ch_references.fasta_fai                ?: Channel.empty() )
    ch_genome_fai_meta              = params.fasta_fai                     ? Channel.fromPath(params.fasta_fai).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                           : ( ch_references.fasta_fai_meta                ?: Channel.empty() )
    ch_mt_shift_fai                 = params.mt_fai_shift                  ? Channel.fromPath(params.mt_fai_shift).collect()
                                                                           : ( ch_references.fasta_fai_mt_shift       ?: Channel.empty() )
    ch_gnomad_af_idx                = params.gnomad_af_idx                 ? Channel.fromPath(params.gnomad_af_idx).collect()
                                                                           : ( ch_references.gnomad_af_idx            ?: Channel.empty() )
    ch_gnomad_af                    = params.gnomad_af                     ? ch_gnomad_af_tab.join(ch_gnomad_af_idx).map {meta, tab, idx -> [tab,idx]}.collect()
                                                                           : Channel.empty()
    ch_gnomad_vcf                   = params.gnomad_vcf                    ? ch_references.gnomad_vcf
                                                                           : Channel.value([])
    ch_known_dbsnp_tbi              = params.known_dbsnp_tbi               ? Channel.fromPath(params.known_dbsnp_tbi).collect()
                                                                           : ( ch_references.known_dbsnp_tbi          ?: Channel.empty() )
    ch_sequence_dictionary          = params.sequence_dictionary           ? Channel.fromPath(params.sequence_dictionary).collect()
                                                                           : ( ch_references.sequence_dict            ?: Channel.empty() )
    ch_sequence_dictionary_mt_shift = params.mt_sequence_dictionary_shift  ? Channel.fromPath(params.mt_sequence_dictionary_shift).collect()
                                                                           : ( ch_references.sequence_dict_mt_shift   ?: Channel.empty() )
    ch_target_bed                   = ch_references.target_bed
    ch_target_intervals             = ch_references.target_intervals
    ch_vcfanno_resources            = params.vcfanno_resources.endsWith('.tar.gz') ? ch_references.vcfanno_resources
                                                                           : Channel.fromPath(params.vcfanno_resources).collect()
    ch_versions                     = ch_versions.mix(ch_references.versions)

    // CREATE CHROMOSOME BED AND INTERVALS
    SCATTER_GENOME (
        ch_sequence_dictionary,
        ch_genome_fai_meta,
        ch_genome_fai_no_meta,
        ch_genome_fasta_no_meta
    )
    .set { ch_scatter }

    ch_scatter_split_intervals  = ch_scatter.split_intervals  ?: Channel.empty()

    // ALIGNING READS, FETCH STATS, AND MERGE.
    ALIGN (
        params.aligner,
        CHECK_INPUT.out.reads,
        ch_genome_fasta_no_meta,
        ch_genome_fai_no_meta,
        ch_bwa_index,
        ch_bwamem2_index,
        ch_known_dbsnp,
        ch_known_dbsnp_tbi
    )
    .set { ch_mapped }
    ch_versions   = ch_versions.mix(ALIGN.out.versions)

    // STEP 1.5: BAM QUALITY CHECK
    QC_BAM (
        ch_mapped.marked_bam,
        ch_mapped.marked_bai,
        ch_genome_fasta_no_meta,
        ch_genome_fai_no_meta,
        ch_bait_intervals,
        ch_target_intervals,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(QC_BAM.out.versions.ifEmpty(null))

    // STEP 1.6: EXPANSIONHUNTER AND STRANGER
    CALL_REPEAT_EXPANSIONS (
        ch_mapped.bam_bai,
        ch_genome_fasta_no_meta,
        ch_variant_catalog
    )
    ch_versions = ch_versions.mix(CALL_REPEAT_EXPANSIONS.out.versions.ifEmpty(null))

    // STEP 2: VARIANT CALLING
    // TODO: There should be a conditional to execute certain variant callers (e.g. sentieon, gatk, deepvariant) defined by the user and we need to think of a default caller.
    CALL_SNV (
        params.variant_caller,
        ch_mapped.bam_bai,
        ch_genome_fasta_no_meta,
        ch_genome_fai_no_meta,
        ch_known_dbsnp,
        ch_known_dbsnp_tbi,
        ch_call_interval,
        ch_ml_model,
        CHECK_INPUT.out.case_info
    )
    ch_versions = ch_versions.mix(CALL_SNV.out.versions)

    CALL_STRUCTURAL_VARIANTS (
        ch_mapped.marked_bam,
        ch_mapped.marked_bai,
        ch_bwa_index,
        ch_genome_fasta_no_meta,
        ch_genome_fasta_meta,
        ch_genome_fai_no_meta,
        CHECK_INPUT.out.case_info,
        ch_target_bed,
        params.cnvpytor_binsizes
    )
    ch_versions = ch_versions.mix(CALL_STRUCTURAL_VARIANTS.out.versions)

    // STEP 2.1: GENS
    if (params.gens_switch) {
        GENS (
            ch_mapped.bam_bai,
            CALL_SNV.out.vcf,
            ch_genome_fasta_meta,
            ch_genome_fai_no_meta,
            file(params.gens_interval_list),
            file(params.gens_pon),
            file(params.gens_gnomad_pos),
            CHECK_INPUT.out.case_info,
            ch_sequence_dictionary
        )
        ch_versions = ch_versions.mix(GENS.out.versions.ifEmpty(null))
    }

    if (params.annotate_sv_switch) {
        ANNOTATE_STRUCTURAL_VARIANTS (
            CALL_STRUCTURAL_VARIANTS.out.vcf,
            params.svdb_query_dbs,
            params.genome,
            params.vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta_no_meta,
            ch_sequence_dictionary
        ).set {ch_sv_annotate}
        ch_versions = ch_versions.mix(ch_sv_annotate.versions)

        ANN_CSQ_PLI_SV (
            ch_sv_annotate.vcf_ann,
            ch_variant_consequences
        )

        RANK_VARIANTS_SV (
            ANN_CSQ_PLI_SV.out.vcf_ann,
            MAKE_PED.out.ped,
            ch_reduced_penetrance,
            ch_score_config_sv
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_SV.out.versions)

        FILTER_VEP_SV(
            RANK_VARIANTS_SV.out.vcf,
            ch_vep_filters
        )

    }

    ANALYSE_MT (
        ch_mapped.bam_bai,
        ch_bwamem2_index,
        ch_genome_fasta_no_meta,
        ch_genome_fasta_meta,
        ch_sequence_dictionary,
        ch_genome_fai_no_meta,
        ch_mt_intervals,
        ch_bwamem2_index_mt_shift,
        ch_mt_fasta_shift_no_meta,
        ch_sequence_dictionary_mt_shift,
        ch_mt_shift_fai,
        ch_mt_intervals_shift,
        ch_mt_backchain_shift,
        params.genome,
        params.vep_cache_version,
        ch_vep_cache,
        CHECK_INPUT.out.case_info
    )
    ch_versions = ch_versions.mix(ANALYSE_MT.out.versions)

    // STEP 3: VARIANT ANNOTATION
    ch_vcf = CALL_SNV.out.vcf.join(CALL_SNV.out.tabix, by: [0])

    if (params.annotate_snv_switch) {
        ANNOTATE_SNVS (
            ch_vcf,
            ch_vcfanno_resources,
            params.vcfanno_toml,
            params.genome,
            params.vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta_no_meta,
            ch_gnomad_af,
            ch_scatter_split_intervals,
            CHECK_INPUT.out.samples
        ).set {ch_snv_annotate}
        ch_versions = ch_versions.mix(ch_snv_annotate.versions)

        ANN_CSQ_PLI_SNV (
            ch_snv_annotate.vcf_ann,
            ch_variant_consequences
        )

        RANK_VARIANTS_SNV (
            ANN_CSQ_PLI_SNV.out.vcf_ann,
            MAKE_PED.out.ped,
            ch_reduced_penetrance,
            ch_score_config_snv
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_SNV.out.versions)

        FILTER_VEP_SNV(
            RANK_VARIANTS_SNV.out.vcf,
            ch_vep_filters
        )

    }

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowRaredisease.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)
    methods_description    = WorkflowRaredisease.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
