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
params.gnomad               = WorkflowMain.getGenomeAttribute(params, 'gnomad')
params.known_dbsnp          = WorkflowMain.getGenomeAttribute(params, 'known_dbsnp')
params.known_indels         = WorkflowMain.getGenomeAttribute(params, 'known_indels')
params.known_mills          = WorkflowMain.getGenomeAttribute(params, 'known_mills')
params.target_bed           = WorkflowMain.getGenomeAttribute(params, 'target_bed')
params.sentieonbwa_index    = WorkflowMain.getGenomeAttribute(params, 'bwa_index')
params.svdb_query_dbs       = WorkflowMain.getGenomeAttribute(params, 'svdb_query_dbs')
params.variant_catalog      = WorkflowMain.getGenomeAttribute(params, 'variant_catalog')
params.vcfanno_resources    = WorkflowMain.getGenomeAttribute(params, 'vcfanno_resources')
params.vcfanno_toml         = WorkflowMain.getGenomeAttribute(params, 'vcfanno_toml')
params.vep_cache            = WorkflowMain.getGenomeAttribute(params, 'vep_cache')

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
