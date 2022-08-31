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

params.fasta                = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fasta_fai            = WorkflowMain.getGenomeAttribute(params, 'fai')
params.bwamem2_index        = WorkflowMain.getGenomeAttribute(params, 'bwamem2')
params.call_interval        = WorkflowMain.getGenomeAttribute(params, 'call_interval')
params.gnomad               = WorkflowMain.getGenomeAttribute(params, 'gnomad')
params.gnomad_af            = WorkflowMain.getGenomeAttribute(params, 'gnomad_af')
params.gnomad_af_tbi        = WorkflowMain.getGenomeAttribute(params, 'gnomad_af_tbi')
params.intervals_mt         = WorkflowMain.getGenomeAttribute(params, 'intervals_mt')
params.known_dbsnp          = WorkflowMain.getGenomeAttribute(params, 'known_dbsnp')
params.known_dbsnp_tbi      = WorkflowMain.getGenomeAttribute(params, 'known_dbsnp_tbi')
params.known_indels         = WorkflowMain.getGenomeAttribute(params, 'known_indels')
params.known_mills          = WorkflowMain.getGenomeAttribute(params, 'known_mills')
params.ml_model             = WorkflowMain.getGenomeAttribute(params, 'ml_model')
params.reduced_penetrance   = WorkflowMain.getGenomeAttribute(params, 'reduced_penetrance')
params.score_config_snv     = WorkflowMain.getGenomeAttribute(params, 'score_config_snv')
params.score_config_sv      = WorkflowMain.getGenomeAttribute(params, 'score_config_sv')
params.target_bed           = WorkflowMain.getGenomeAttribute(params, 'target_bed')
params.sentieonbwa_index    = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.svdb_query_dbs       = WorkflowMain.getGenomeAttribute(params, 'svdb_query_dbs')
params.variant_catalog      = WorkflowMain.getGenomeAttribute(params, 'variant_catalog')
params.vcfanno_resources    = WorkflowMain.getGenomeAttribute(params, 'vcfanno_resources')
params.vcfanno_toml         = WorkflowMain.getGenomeAttribute(params, 'vcfanno_toml')
params.vep_cache            = WorkflowMain.getGenomeAttribute(params, 'vep_cache')
params.gens_interval_list   = WorkflowMain.getGenomeAttribute(params, 'gens_interval_list')
params.gens_pon             = WorkflowMain.getGenomeAttribute(params, 'gens_pon')
params.gens_gnomad_pos      = WorkflowMain.getGenomeAttribute(params, 'gens_gnomad_pos')

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
