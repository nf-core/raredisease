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
params.fasta_fai                      = WorkflowMain.getGenomeAttribute(params, 'fai')
params.bwa                            = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.bwamem2                        = WorkflowMain.getGenomeAttribute(params, 'bwamem2')
params.call_interval                  = WorkflowMain.getGenomeAttribute(params, 'call_interval')
params.gnomad_af                      = WorkflowMain.getGenomeAttribute(params, 'gnomad_af')
params.gnomad_af_idx                  = WorkflowMain.getGenomeAttribute(params, 'gnomad_af_idx')
params.intervals_wgs                  = WorkflowMain.getGenomeAttribute(params, 'intervals_wgs')
params.intervals_y                    = WorkflowMain.getGenomeAttribute(params, 'intervals_y')
params.known_dbsnp                    = WorkflowMain.getGenomeAttribute(params, 'known_dbsnp')
params.known_dbsnp_tbi                = WorkflowMain.getGenomeAttribute(params, 'known_dbsnp_tbi')
params.known_indels                   = WorkflowMain.getGenomeAttribute(params, 'known_indels')
params.known_mills                    = WorkflowMain.getGenomeAttribute(params, 'known_mills')
params.ml_model                       = WorkflowMain.getGenomeAttribute(params, 'ml_model')
params.mt_backchain_shift             = WorkflowMain.getGenomeAttribute(params, 'mt_backchain_shift')
params.mt_bwa_index_shift             = WorkflowMain.getGenomeAttribute(params, 'mt_bwa_index_shift')
params.mt_bwamem2_index_shift         = WorkflowMain.getGenomeAttribute(params, 'mt_bwamem2_index_shift')
params.mt_fasta_shift                 = WorkflowMain.getGenomeAttribute(params, 'mt_fasta_shift')
params.mt_fai_shift                   = WorkflowMain.getGenomeAttribute(params, 'mt_fai_shift')
params.mt_intervals                   = WorkflowMain.getGenomeAttribute(params, 'mt_intervals')
params.mt_intervals_shift             = WorkflowMain.getGenomeAttribute(params, 'mt_intervals_shift')
params.mt_sequence_dictionary_shift   = WorkflowMain.getGenomeAttribute(params, 'mt_sequence_dictionary_shift')
params.reduced_penetrance             = WorkflowMain.getGenomeAttribute(params, 'reduced_penetrance')
params.sequence_dictionary            = WorkflowMain.getGenomeAttribute(params, 'sequence_dictionary')
params.score_config_snv               = WorkflowMain.getGenomeAttribute(params, 'score_config_snv')
params.score_config_sv                = WorkflowMain.getGenomeAttribute(params, 'score_config_sv')
params.target_bed                     = WorkflowMain.getGenomeAttribute(params, 'target_bed')
params.svdb_query_dbs                 = WorkflowMain.getGenomeAttribute(params, 'svdb_query_dbs')
params.variant_catalog                = WorkflowMain.getGenomeAttribute(params, 'variant_catalog')
params.vep_filters                    = WorkflowMain.getGenomeAttribute(params, 'vep_filters')
params.vcfanno_resources              = WorkflowMain.getGenomeAttribute(params, 'vcfanno_resources')
params.vcfanno_toml                   = WorkflowMain.getGenomeAttribute(params, 'vcfanno_toml')
params.vcfanno_lua                    = WorkflowMain.getGenomeAttribute(params, 'vcfanno_lua')
params.vep_cache                      = WorkflowMain.getGenomeAttribute(params, 'vep_cache')
params.vep_cache_version              = WorkflowMain.getGenomeAttribute(params, 'vep_cache_version')
params.gens_interval_list             = WorkflowMain.getGenomeAttribute(params, 'gens_interval_list')
params.gens_pon                       = WorkflowMain.getGenomeAttribute(params, 'gens_pon')
params.gens_gnomad_pos                = WorkflowMain.getGenomeAttribute(params, 'gens_gnomad_pos')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

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
