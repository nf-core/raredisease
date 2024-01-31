#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/raredisease
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/raredisease
    Website: https://nf-co.re/raredisease
    Slack  : https://nfcore.slack.com/channels/raredisease
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.fasta                          = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fai                            = WorkflowMain.getGenomeAttribute(params, 'fai')
params.bwa                            = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.bwamem2                        = WorkflowMain.getGenomeAttribute(params, 'bwamem2')
params.call_interval                  = WorkflowMain.getGenomeAttribute(params, 'call_interval')
params.cadd_resources                 = WorkflowMain.getGenomeAttribute(params, 'cadd_resources')
params.gcnvcaller_model               = WorkflowMain.getGenomeAttribute(params, 'gcnvcaller_model')
params.gens_interval_list             = WorkflowMain.getGenomeAttribute(params, 'gens_interval_list')
params.gens_pon                       = WorkflowMain.getGenomeAttribute(params, 'gens_pon')
params.gens_gnomad_pos                = WorkflowMain.getGenomeAttribute(params, 'gens_gnomad_pos')
params.gnomad_af                      = WorkflowMain.getGenomeAttribute(params, 'gnomad_af')
params.gnomad_af_idx                  = WorkflowMain.getGenomeAttribute(params, 'gnomad_af_idx')
params.intervals_wgs                  = WorkflowMain.getGenomeAttribute(params, 'intervals_wgs')
params.intervals_y                    = WorkflowMain.getGenomeAttribute(params, 'intervals_y')
params.known_dbsnp                    = WorkflowMain.getGenomeAttribute(params, 'known_dbsnp')
params.known_dbsnp_tbi                = WorkflowMain.getGenomeAttribute(params, 'known_dbsnp_tbi')
params.mobile_element_references      = WorkflowMain.getGenomeAttribute(params, 'mobile_element_references')
params.ml_model                       = WorkflowMain.getGenomeAttribute(params, 'ml_model')
params.mt_fasta                       = WorkflowMain.getGenomeAttribute(params, 'mt_fasta')
params.ploidy_model                   = WorkflowMain.getGenomeAttribute(params, 'ploidy_model')
params.reduced_penetrance             = WorkflowMain.getGenomeAttribute(params, 'reduced_penetrance')
params.readcount_intervals            = WorkflowMain.getGenomeAttribute(params, 'readcount_intervals')
params.rtg_truthvcfs                  = WorkflowMain.getGenomeAttribute(params, 'rtg_truthvcfs')
params.sequence_dictionary            = WorkflowMain.getGenomeAttribute(params, 'sequence_dictionary')
params.score_config_mt                = WorkflowMain.getGenomeAttribute(params, 'score_config_mt')
params.score_config_snv               = WorkflowMain.getGenomeAttribute(params, 'score_config_snv')
params.score_config_sv                = WorkflowMain.getGenomeAttribute(params, 'score_config_sv')
params.sdf                            = WorkflowMain.getGenomeAttribute(params, 'sdf')
params.svdb_query_dbs                 = WorkflowMain.getGenomeAttribute(params, 'svdb_query_dbs')
params.target_bed                     = WorkflowMain.getGenomeAttribute(params, 'target_bed')
params.variant_catalog                = WorkflowMain.getGenomeAttribute(params, 'variant_catalog')
params.vep_filters                    = WorkflowMain.getGenomeAttribute(params, 'vep_filters')
params.vcf2cytosure_blacklist         = WorkflowMain.getGenomeAttribute(params, 'vcf2cytosure_blacklist')
params.vcfanno_resources              = WorkflowMain.getGenomeAttribute(params, 'vcfanno_resources')
params.vcfanno_toml                   = WorkflowMain.getGenomeAttribute(params, 'vcfanno_toml')
params.vcfanno_lua                    = WorkflowMain.getGenomeAttribute(params, 'vcfanno_lua')
params.vep_cache                      = WorkflowMain.getGenomeAttribute(params, 'vep_cache')
params.vep_plugin_files               = WorkflowMain.getGenomeAttribute(params, 'vep_plugin_files')
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

// Print help message if needed
if (params.help) {
    def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
    def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
    def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
    System.exit(0)
}

// Validate input parameters
if (params.validate_params) {
    validateParameters()
}

WorkflowMain.initialise(workflow, params, log, args)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RAREDISEASE } from './workflows/raredisease'

//
// WORKFLOW: Run main nf-core/raredisease analysis pipeline
//
workflow NFCORE_RAREDISEASE {
    RAREDISEASE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_RAREDISEASE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
