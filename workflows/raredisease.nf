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
    params.bwamem2_index,
    params.fasta,
    params.fasta_fai,
    params.gnomad,
    params.input,
    params.intervals_mt,
    params.multiqc_config,
    params.reduced_penetrance,
    params.score_config_snv,
    params.score_config_sv,
    params.sentieonbwa_index,
    params.svdb_query_dbs,
    params.vcfanno_resources
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
ch_ml_model             = params.ml_model            ? file(params.ml_model)           : []
ch_call_interval        = params.call_interval       ? file(params.call_interval)      : []
ch_reduced_penetrance   = params.reduced_penetrance  ? file(params.reduced_penetrance) : []
ch_score_config_snv     = params.score_config_snv    ? file(params.score_config_snv)   : []
ch_score_config_sv      = params.score_config_sv     ? file(params.score_config_sv)    : []
ch_variant_consequences = file("$projectDir/assets/variant_consequences_v1.txt", checkIfExists: true)
ch_vep_cache            = params.vep_cache           ? file(params.vep_cache)          : []
ch_vep_filters          = params.vep_filters         ? file(params.vep_filters)        : []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: local modules
//

include { MAKE_PED                     } from '../modules/local/create_pedfile'
include { FILTER_VEP as FILTER_VEP_SNV } from '../modules/local/filter_vep'
include { FILTER_VEP as FILTER_VEP_SV  } from '../modules/local/filter_vep'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { CHECK_INPUT                  } from '../subworkflows/local/check_input'
include { PREPARE_REFERENCES           } from '../subworkflows/local/prepare_references'
include { ANNOTATE_SNVS                } from '../subworkflows/local/annotate_snvs'
include { ANNOTATE_STRUCTURAL_VARIANTS } from '../subworkflows/local/annotate_structural_variants'
include { GENS                         } from '../subworkflows/local/gens'
include { ALIGN                        } from '../subworkflows/local/align'
include { CALL_SNV                     } from '../subworkflows/local/call_snv'
include { ANALYSE_MT                   } from '../subworkflows/local/analyse_MT'
include { ANNOTATE_CSQ as ANN_CSQ_SNV  } from '../subworkflows/local/annotate_consequence'
include { ANNOTATE_CSQ as ANN_CSQ_SV   } from '../subworkflows/local/annotate_consequence'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

//
// SUBWORKFLOW: Consists entirely of nf-core/modules
//

include { CALL_REPEAT_EXPANSIONS             } from '../subworkflows/nf-core/call_repeat_expansions'
include { QC_BAM                             } from '../subworkflows/nf-core/qc_bam'
include { CALL_STRUCTURAL_VARIANTS           } from '../subworkflows/nf-core/call_structural_variants'
include { RANK_VARIANTS as RANK_VARIANTS_SNV } from '../subworkflows/nf-core/genmod'
include { RANK_VARIANTS as RANK_VARIANTS_SV  } from '../subworkflows/nf-core/genmod'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RAREDISEASE {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    CHECK_INPUT (
        ch_input
    )
    ch_versions = ch_versions.mix(CHECK_INPUT.out.versions)

    MAKE_PED (CHECK_INPUT.out.samples.toList())

    // STEP 0: QUALITY CHECK.
    FASTQC (
        CHECK_INPUT.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // STEP 0: PREPARE GENOME REFERENCES AND INDICES.
    PREPARE_REFERENCES (
        params.aligner,
        params.bwamem2_index,
        params.fasta,
        params.fasta_fai,
        params.gnomad,
        params.gnomad_af,
        params.gnomad_af_tbi,
        params.known_dbsnp,
        params.known_dbsnp_tbi,
        params.sentieonbwa_index,
        params.target_bed,
        params.variant_catalog,
        params.vcfanno_resources
    )
    .set { ch_references }
    ch_versions = ch_versions.mix(ch_references.versions)

    // STEP 1: ALIGNING READS, FETCH STATS, AND MERGE.
    ALIGN (
        params.aligner,
        CHECK_INPUT.out.reads,
        ch_references.genome_fasta,
        ch_references.genome_fai,
        ch_references.aligner_index,
        ch_references.known_dbsnp,
        ch_references.known_dbsnp_tbi
    )
    .set { ch_mapped }
    ch_versions   = ch_versions.mix(ALIGN.out.versions)

    // STEP 1.5: BAM QUALITY CHECK
    QC_BAM (
        ch_mapped.marked_bam,
        ch_mapped.marked_bai,
        ch_references.genome_fasta,
        ch_references.genome_fai,
        ch_references.bait_intervals,
        ch_references.target_intervals,
        ch_references.chrom_sizes
    )
    ch_versions = ch_versions.mix(QC_BAM.out.versions.ifEmpty(null))

    // STEP 1.6: EXPANSIONHUNTER AND STRANGER
    CALL_REPEAT_EXPANSIONS (
        ch_mapped.bam_bai,
        ch_references.genome_fasta,
        ch_references.variant_catalog
    )
    ch_versions = ch_versions.mix(CALL_REPEAT_EXPANSIONS.out.versions.ifEmpty(null))

    // STEP 2: VARIANT CALLING
    // TODO: There should be a conditional to execute certain variant callers (e.g. sentieon, gatk, deepvariant) defined by the user and we need to think of a default caller.
    CALL_SNV (
        params.variant_caller,
        ch_mapped.bam_bai,
        ch_references.genome_fasta,
        ch_references.genome_fai,
        ch_references.known_dbsnp,
        ch_references.known_dbsnp_tbi,
        ch_call_interval,
        ch_ml_model,
        CHECK_INPUT.out.case_info
    )
    ch_versions = ch_versions.mix(CALL_SNV.out.versions)

    CALL_STRUCTURAL_VARIANTS (
        ch_mapped.marked_bam,
        ch_mapped.marked_bai,
        ch_references.genome_fasta,
        ch_references.genome_fai,
        CHECK_INPUT.out.case_info,
        ch_references.target_bed,
        params.cnvpytor_binsizes
    )
    ch_versions = ch_versions.mix(CALL_STRUCTURAL_VARIANTS.out.versions)

    // STEP 2.1: GENS
    if (params.gens_switch) {
        GENS (
            ch_mapped.bam_bai,
            CALL_SNV.out.vcf,
            ch_references.genome_fasta,
            ch_references.genome_fai,
            file(params.gens_interval_list),
            file(params.gens_pon),
            file(params.gens_gnomad_pos),
            CHECK_INPUT.out.case_info,
            ch_references.sequence_dict
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
            ch_references.genome_fasta,
            ch_references.sequence_dict
        ).set {ch_sv_annotate}
        ch_versions = ch_versions.mix(ch_sv_annotate.versions)

        RANK_VARIANTS_SV (
            ch_sv_annotate.vcf_ann,
            MAKE_PED.out.ped,
            ch_reduced_penetrance,
            ch_score_config_sv
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_SV.out.versions)

        FILTER_VEP_SV(
            RANK_VARIANTS_SV.out.vcf,
            ch_vep_filters
        )

        ANN_CSQ_SV (
            FILTER_VEP_SV.out.vcf,
            RANK_VARIANTS_SV.out.vcf,
            ch_variant_consequences
        )
    }


    // STEP 2.1: ANALYSE MT
    ch_intervals_mt = Channel.fromPath(params.intervals_mt)
    ANALYSE_MT (
        ch_mapped.bam_bai,
        ch_references.aligner_index,
        ch_references.genome_fasta,
        ch_references.sequence_dict,
        ch_references.genome_fai,
        ch_intervals_mt
    )
    ch_versions = ch_versions.mix(ANALYSE_MT.out.versions)

    // STEP 3: VARIANT ANNOTATION
    ch_vcf = CALL_SNV.out.vcf.join(CALL_SNV.out.tabix, by: [0])

    if (params.annotate_snv_switch) {
        ANNOTATE_SNVS (
            ch_vcf,
            ch_references.vcfanno_resources,
            params.vcfanno_toml,
            params.genome,
            params.vep_cache_version,
            ch_vep_cache,
            ch_references.genome_fasta,
            ch_references.gnomad_af,
            CHECK_INPUT.out.samples
        ).set {ch_snv_annotate}
        ch_versions = ch_versions.mix(ch_snv_annotate.versions)

        RANK_VARIANTS_SNV (
            ch_snv_annotate.vcf_ann,
            MAKE_PED.out.ped,
            ch_reduced_penetrance,
            ch_score_config_snv
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_SNV.out.versions)

        FILTER_VEP_SNV(
            RANK_VARIANTS_SNV.out.vcf,
            ch_vep_filters
        )

        ANN_CSQ_SNV (
            FILTER_VEP_SNV.out.vcf,
            RANK_VARIANTS_SNV.out.vcf,
            ch_variant_consequences
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

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        [],
        [],
        []
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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
