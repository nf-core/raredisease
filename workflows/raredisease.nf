/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.bwa,
    params.bwamem2,
    params.call_interval,
    params.cadd_resources,
    params.fasta,
    params.fai,
    params.gens_gnomad_pos,
    params.gens_interval_list,
    params.gens_pon,
    params.gnomad_af,
    params.gnomad_af_idx,
    params.input,
    params.intervals_wgs,
    params.intervals_y,
    params.known_dbsnp,
    params.known_dbsnp_tbi,
    params.known_indels,
    params.known_mills,
    params.ml_model,
    params.mt_fasta,
    params.multiqc_config,
    params.reduced_penetrance,
    params.score_config_snv,
    params.score_config_sv,
    params.sequence_dictionary,
    params.target_bed,
    params.svdb_query_dbs,
    params.variant_catalog,
    params.vep_filters,
    params.vcfanno_lua,
    params.vcfanno_resources,
    params.vcfanno_toml,
    params.vep_cache
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowRaredisease.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHECK MANDATORY PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def mandatoryParams = [
    "aligner",
    "analysis_type",
    "fasta",
    "input",
    "intervals_wgs",
    "intervals_y",
    "platform",
    "variant_catalog",
    "variant_caller"
]

if (!params.skip_snv_annotation) {
    mandatoryParams += ["genome", "vcfanno_resources", "vcfanno_toml", "vep_cache", "vep_cache_version",
    "gnomad_af", "score_config_snv"]
}

if (!params.skip_sv_annotation) {
    mandatoryParams += ["genome", "svdb_query_dbs", "vep_cache", "vep_cache_version", "score_config_sv"]
}

if (!params.skip_mt_analysis) {
    mandatoryParams += ["genome", "mito_name", "vcfanno_resources", "vcfanno_toml", "vep_cache_version", "vep_cache"]
}

if (params.analysis_type.equals("wes")) {
    mandatoryParams += ["target_bed"]
}

if (params.variant_caller.equals("sentieon")) {
    mandatoryParams += ["ml_model"]
}

def missingParamsCount = 0
for (param in mandatoryParams.unique()) {
    if (params[param] == null) {
        println("params." + param + " not set.")
        missingParamsCount += 1
    }
}

if (missingParamsCount>0) {
    error("\nSet missing parameters and restart the run.")
}
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

//
// MODULE: Installed directly from nf-core/modules
//

include { BCFTOOLS_CONCAT                       } from '../modules/nf-core/bcftools/concat/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS           } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTQC                                } from '../modules/nf-core/fastqc/main'
include { GATK4_SELECTVARIANTS                  } from '../modules/nf-core/gatk4/selectvariants/main'
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'
include { SMNCOPYNUMBERCALLER                   } from '../modules/nf-core/smncopynumbercaller/main'

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
include { PEDDY_CHECK                           } from '../subworkflows/local/peddy_check'


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
    ch_input = Channel.fromPath(params.input)
    CHECK_INPUT (ch_input)
    ch_versions = ch_versions.mix(CHECK_INPUT.out.versions)

    // Initialize file channels for PREPARE_REFERENCES subworkflow
    ch_genome_fasta             = Channel.fromPath(params.fasta).map { it -> [[id:it[0].simpleName], it] }.collect()
    ch_genome_fai               = params.fai            ? Channel.fromPath(params.fai).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                        : Channel.empty()
    ch_gnomad_af_tab            = params.gnomad_af      ? Channel.fromPath(params.gnomad_af).map{ it -> [[id:it[0].simpleName], it] }.collect()
                                                        : Channel.value([[],[]])
    ch_dbsnp                    = params.known_dbsnp    ? Channel.fromPath(params.known_dbsnp).map{ it -> [[id:it[0].simpleName], it] }.collect()
                                                        : Channel.value([[],[]])
    ch_mt_fasta                 = params.mt_fasta       ? Channel.fromPath(params.mt_fasta).map { it -> [[id:it[0].simpleName], it] }.collect()
                                                        : Channel.empty()
    ch_target_bed_unprocessed   = params.target_bed     ? Channel.fromPath(params.target_bed).map{ it -> [[id:it[0].simpleName], it] }.collect()
                                                        : Channel.value([[],[]])
    ch_vep_cache_unprocessed    = params.vep_cache      ? Channel.fromPath(params.vep_cache).map { it -> [[id:'vep_cache'], it] }.collect()
                                                        : Channel.value([[],[]])

    // Prepare references and indices.
    PREPARE_REFERENCES (
        ch_genome_fasta,
        ch_genome_fai,
        ch_mt_fasta,
        ch_gnomad_af_tab,
        ch_dbsnp,
        ch_target_bed_unprocessed,
        ch_vep_cache_unprocessed
    )
    .set { ch_references }


    // Gather built indices or get them from the params
    ch_bait_intervals           = ch_references.bait_intervals
    ch_cadd_header              = Channel.fromPath("$projectDir/assets/cadd_to_vcf_header_-1.0-.txt", checkIfExists: true).collect()
    ch_cadd_resources           = params.cadd_resources                    ? Channel.fromPath(params.cadd_resources).collect()
                                                                           : Channel.value([])
    ch_call_interval            = params.call_interval                     ? Channel.fromPath(params.call_interval).collect()
                                                                           : Channel.value([])
    ch_dbsnp_tbi                = params.known_dbsnp_tbi                   ? Channel.fromPath(params.known_dbsnp_tbi).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                           : ch_references.known_dbsnp_tbi.ifEmpty([[],[]])
    ch_genome_bwaindex          = params.bwa                               ? Channel.fromPath(params.bwa).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                           : ch_references.genome_bwa_index
    ch_genome_bwamem2index      = params.bwamem2                           ? Channel.fromPath(params.bwamem2).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                           : ch_references.genome_bwamem2_index
    ch_genome_chrsizes          = ch_references.genome_chrom_sizes
    ch_genome_fai               = ch_references.genome_fai
    ch_genome_dictionary        = params.sequence_dictionary               ? Channel.fromPath(params.sequence_dictionary).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                           : ch_references.genome_dict
    ch_gnomad_afidx             = params.gnomad_af_idx                     ? Channel.fromPath(params.gnomad_af_idx).collect()
                                                                           : ch_references.gnomad_af_idx
    ch_gnomad_af                = params.gnomad_af                         ? ch_gnomad_af_tab.join(ch_gnomad_afidx).map {meta, tab, idx -> [tab,idx]}.collect()
                                                                           : Channel.empty()
    ch_intervals_wgs            = params.intervals_wgs                     ? Channel.fromPath(params.intervals_wgs).collect()
                                                                           : Channel.empty()
    ch_intervals_y              = params.intervals_y                       ? Channel.fromPath(params.intervals_y).collect()
                                                                           : Channel.empty()
    ch_ml_model                 = params.variant_caller.equals("sentieon") ? Channel.fromPath(params.ml_model).collect()
                                                                           : Channel.value([])
    ch_mt_intervals             = ch_references.mt_intervals
    ch_mtshift_backchain        = ch_references.mtshift_backchain
    ch_mtshift_bwaindex         = ch_references.mtshift_bwa_index
    ch_mtshift_bwamem2index     = ch_references.mtshift_bwamem2_index
    ch_mtshift_dictionary       = ch_references.mtshift_dict
    ch_mtshift_fai              = ch_references.mtshift_fai
    ch_mtshift_fasta            = ch_references.mtshift_fasta
    ch_mtshift_intervals        = ch_references.mtshift_intervals
    ch_reduced_penetrance       = params.reduced_penetrance                ? Channel.fromPath(params.reduced_penetrance).collect()
                                                                           : Channel.value([])
    ch_score_config_snv         = params.score_config_snv                  ? Channel.fromPath(params.score_config_snv).collect()
                                                                           : Channel.value([])
    ch_score_config_sv          = params.score_config_sv                   ? Channel.fromPath(params.score_config_sv).collect()
                                                                           : Channel.value([])
    ch_target_bed               = ch_references.target_bed
    ch_target_intervals         = ch_references.target_intervals
    ch_variant_catalog          = params.variant_catalog                   ? Channel.fromPath(params.variant_catalog).map { it -> [[id:it[0].simpleName],it]}.collect()
                                                                           : Channel.value([[],[]])
    ch_variant_consequences     = Channel.fromPath("$projectDir/assets/variant_consequences_v1.txt", checkIfExists: true).collect()
    ch_vcfanno_resources        = params.vcfanno_resources                 ? Channel.fromPath(params.vcfanno_resources).splitText().map{it -> it.trim()}.collect()
                                                                           : Channel.value([])
    ch_vcfanno_lua              = params.vcfanno_lua                       ? Channel.fromPath(params.vcfanno_lua).collect()
                                                                           : Channel.value([])
    ch_vcfanno_toml             = params.vcfanno_toml                      ? Channel.fromPath(params.vcfanno_toml).collect()
                                                                           : Channel.value([])
    ch_vep_cache                = ( params.vep_cache && params.vep_cache.endsWith("tar.gz") )  ? ch_references.vep_resources
                                                                           : ( params.vep_cache    ? Channel.fromPath(params.vep_cache).collect() : Channel.value([]) )
    ch_vep_filters              = params.vep_filters                       ? Channel.fromPath(params.vep_filters).collect()
                                                                           : Channel.value([])
    ch_versions                 = ch_versions.mix(ch_references.versions)


    // Generate pedigree file
    ch_pedfile = CHECK_INPUT.out.samples.toList().map { makePed(it) }

    // Input QC
    FASTQC (CHECK_INPUT.out.reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // CREATE CHROMOSOME BED AND INTERVALS
    SCATTER_GENOME (
        ch_genome_dictionary,
        ch_genome_fai,
        ch_genome_fasta
    )
    .set { ch_scatter }

    ch_scatter_split_intervals  = ch_scatter.split_intervals  ?: Channel.empty()

    // ALIGNING READS, FETCH STATS, AND MERGE.
    ALIGN (
        CHECK_INPUT.out.reads,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_bwaindex,
        ch_genome_bwamem2index,
        ch_dbsnp,
        ch_dbsnp_tbi,
        params.platform
    )
    .set { ch_mapped }
    ch_versions   = ch_versions.mix(ALIGN.out.versions)

    // BAM QUALITY CHECK
    QC_BAM (
        ch_mapped.marked_bam,
        ch_mapped.marked_bai,
        ch_mapped.bam_bai,
        ch_genome_fasta,
        ch_genome_fai,
        ch_bait_intervals,
        ch_target_intervals,
        ch_genome_chrsizes,
        ch_intervals_wgs,
        ch_intervals_y
    )
    ch_versions = ch_versions.mix(QC_BAM.out.versions)

    // EXPANSIONHUNTER AND STRANGER
    CALL_REPEAT_EXPANSIONS (
        ch_mapped.bam_bai,
        ch_variant_catalog,
        CHECK_INPUT.out.case_info,
        ch_genome_fasta,
        ch_genome_fai
    )
    ch_versions = ch_versions.mix(CALL_REPEAT_EXPANSIONS.out.versions)

    // STEP 1.7: SMNCOPYNUMBERCALLER
    ch_mapped.bam_bai
        .collect{it[1]}
        .toList()
        .set { ch_bam_list }

    ch_mapped.bam_bai
        .collect{it[2]}
        .toList()
        .set { ch_bai_list }

    CHECK_INPUT.out.case_info
        .combine(ch_bam_list)
        .combine(ch_bai_list)
        .set { ch_bams_bais }

    SMNCOPYNUMBERCALLER (
        ch_bams_bais
    )
    ch_versions = ch_versions.mix(SMNCOPYNUMBERCALLER.out.versions)

    // STEP 2: VARIANT CALLING
    CALL_SNV (
        ch_mapped.bam_bai,
        ch_genome_fasta,
        ch_genome_fai,
        ch_dbsnp,
        ch_dbsnp_tbi,
        ch_call_interval,
        ch_ml_model,
        CHECK_INPUT.out.case_info
    )
    ch_versions = ch_versions.mix(CALL_SNV.out.versions)

    CALL_STRUCTURAL_VARIANTS (
        ch_mapped.marked_bam,
        ch_mapped.marked_bai,
        ch_mapped.bam_bai,
        ch_genome_bwaindex,
        ch_genome_fasta,
        ch_genome_fai,
        CHECK_INPUT.out.case_info,
        ch_target_bed
    )
    ch_versions = ch_versions.mix(CALL_STRUCTURAL_VARIANTS.out.versions)

    // ped correspondence, sex check, ancestry check
    PEDDY_CHECK (
        CALL_SNV.out.vcf.join(CALL_SNV.out.tabix, failOnMismatch:true, failOnDuplicate:true),
        ch_pedfile
    )
    ch_versions = ch_versions.mix(PEDDY_CHECK.out.versions)

    // GENS
    if (params.gens_switch) {
        GENS (
            ch_mapped.bam_bai,
            CALL_SNV.out.vcf,
            ch_genome_fasta,
            ch_genome_fai,
            file(params.gens_interval_list),
            file(params.gens_pon),
            file(params.gens_gnomad_pos),
            CHECK_INPUT.out.case_info,
            ch_genome_dictionary
        )
        ch_versions = ch_versions.mix(GENS.out.versions)
    }

    if (!params.skip_sv_annotation) {
        ANNOTATE_STRUCTURAL_VARIANTS (
            CALL_STRUCTURAL_VARIANTS.out.vcf,
            params.svdb_query_dbs,
            params.genome,
            params.vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta,
            ch_genome_dictionary
        ).set {ch_sv_annotate}
        ch_versions = ch_versions.mix(ch_sv_annotate.versions)

        ANN_CSQ_PLI_SV (
            ch_sv_annotate.vcf_ann,
            ch_variant_consequences
        )
        ch_versions = ch_versions.mix(ANN_CSQ_PLI_SV.out.versions)

        RANK_VARIANTS_SV (
            ANN_CSQ_PLI_SV.out.vcf_ann,
            ch_pedfile,
            ch_reduced_penetrance,
            ch_score_config_sv
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_SV.out.versions)

        FILTER_VEP_SV(
            RANK_VARIANTS_SV.out.vcf,
            ch_vep_filters
        )
        ch_versions = ch_versions.mix(FILTER_VEP_SV.out.versions)

    }

    if (!params.skip_mt_analysis) {
        ANALYSE_MT (
            ch_mapped.bam_bai,
            ch_cadd_header,
            ch_cadd_resources,
            ch_genome_bwaindex,
            ch_genome_bwamem2index,
            ch_genome_fasta,
            ch_genome_fai,
            ch_genome_dictionary,
            ch_mt_intervals,
            ch_mtshift_bwaindex,
            ch_mtshift_bwamem2index,
            ch_mtshift_fasta,
            ch_mtshift_dictionary,
            ch_mtshift_fai,
            ch_mtshift_intervals,
            ch_mtshift_backchain,
            ch_vcfanno_resources,
            ch_vcfanno_toml,
            params.genome,
            params.vep_cache_version,
            ch_vep_cache,
            CHECK_INPUT.out.case_info
        )

        ch_versions = ch_versions.mix(ANALYSE_MT.out.versions)

    }

    // VARIANT ANNOTATION

    if (!params.skip_snv_annotation) {

        ch_vcf = CALL_SNV.out.vcf.join(CALL_SNV.out.tabix, failOnMismatch:true, failOnDuplicate:true)

        if (!params.skip_mt_analysis) {
            ch_vcf
                .map { meta, vcf, tbi -> return [meta, vcf, tbi, []]}
                .set { ch_selvar_in }

            GATK4_SELECTVARIANTS(ch_selvar_in) // remove mitochondrial variants

            ch_vcf = GATK4_SELECTVARIANTS.out.vcf.join(GATK4_SELECTVARIANTS.out.tbi, failOnMismatch:true, failOnDuplicate:true)
            ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions)
        }

        ANNOTATE_SNVS (
            ch_vcf,
            params.analysis_type,
            ch_cadd_header,
            ch_cadd_resources,
            ch_vcfanno_resources,
            ch_vcfanno_lua,
            ch_vcfanno_toml,
            params.genome,
            params.vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta,
            ch_gnomad_af,
            ch_scatter_split_intervals
        ).set {ch_snv_annotate}
        ch_versions = ch_versions.mix(ch_snv_annotate.versions)

        ch_snv_annotate = ANNOTATE_SNVS.out.vcf_ann

        if (!params.skip_mt_analysis) {

            ANNOTATE_SNVS.out.vcf_ann
                .concat(ANALYSE_MT.out.vcf)
                .groupTuple()
                .set { ch_merged_vcf }

            ANNOTATE_SNVS.out.tbi
                .concat(ANALYSE_MT.out.tbi)
                .groupTuple()
                .set { ch_merged_tbi }

            ch_merged_vcf.join(ch_merged_tbi, failOnMismatch:true, failOnDuplicate:true).set {ch_concat_in}

            BCFTOOLS_CONCAT (ch_concat_in)
            ch_snv_annotate = BCFTOOLS_CONCAT.out.vcf
            ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
        }

        ANN_CSQ_PLI_SNV (
            ch_snv_annotate,
            ch_variant_consequences
        )
        ch_versions = ch_versions.mix(ANN_CSQ_PLI_SNV.out.versions)

        RANK_VARIANTS_SNV (
            ANN_CSQ_PLI_SNV.out.vcf_ann,
            ch_pedfile,
            ch_reduced_penetrance,
            ch_score_config_snv
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_SNV.out.versions)

        FILTER_VEP_SNV(
            RANK_VARIANTS_SNV.out.vcf,
            ch_vep_filters
        )
        ch_versions = ch_versions.mix(FILTER_VEP_SNV.out.versions)

    }

    //
    // MODULE: Pipeline reporting
    //

    // The template v2.7.1 template update introduced: ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml')
    // This caused the pipeline to stall
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowRaredisease.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowRaredisease.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.multiple_metrics.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.hs_metrics.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.qualimap_results.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.global_dist.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.cov.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PEDDY_CHECK.out.ped.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PEDDY_CHECK.out.csv.map{it[1]}.collect().ifEmpty([]))


    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
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
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def makePed(samples) {

    def case_name  = samples[0].case_id
    def outfile  = file("${params.outdir}/pipeline_info/${case_name}" + '.ped')
    outfile.text = ['#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype'].join('\t')
    def samples_list = []
    for(int i = 0; i<samples.size(); i++) {
        sample_tokenized   =  samples[i].id.tokenize("_")
        sample_tokenized.removeLast()
        sample_name        =  sample_tokenized.join("_")
        if (!samples_list.contains(sample_name)) {
            outfile.append('\n' + [samples[i].case_id, sample_name, samples[i].paternal, samples[i].maternal, samples[i].sex, samples[i].phenotype].join('\t'));
            samples_list.add(sample_name)
        }
    }
    return outfile
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
