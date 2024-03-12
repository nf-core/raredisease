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

params.fasta                           = getGenomeAttribute('fasta')
params.fai                             = getGenomeAttribute('fai')
params.bwa                             = getGenomeAttribute('bwa')
params.bwamem2                         = getGenomeAttribute('bwamem2')
params.call_interval                   = getGenomeAttribute('call_interval')
params.cadd_resources                  = getGenomeAttribute('cadd_resources')
params.gcnvcaller_model                = getGenomeAttribute('gcnvcaller_model')
params.gens_interval_list              = getGenomeAttribute('gens_interval_list')
params.gens_pon_female                 = getGenomeAttribute('gens_pon_female')
params.gens_pon_male                   = getGenomeAttribute('gens_pon_male')
params.gens_gnomad_pos                 = getGenomeAttribute('gens_gnomad_pos')
params.gnomad_af                       = getGenomeAttribute('gnomad_af')
params.gnomad_af_idx                   = getGenomeAttribute('gnomad_af_idx')
params.intervals_wgs                   = getGenomeAttribute('intervals_wgs')
params.intervals_y                     = getGenomeAttribute('intervals_y')
params.known_dbsnp                     = getGenomeAttribute('known_dbsnp')
params.known_dbsnp_tbi                 = getGenomeAttribute('known_dbsnp_tbi')
params.mobile_element_references       = getGenomeAttribute('mobile_element_references')
params.mobile_element_svdb_annotations = getGenomeAttribute('mobile_element_svdb_annotations')
params.ml_model                        = getGenomeAttribute('ml_model')
params.mt_fasta                        = getGenomeAttribute('mt_fasta')
params.ploidy_model                    = getGenomeAttribute('ploidy_model')
params.reduced_penetrance              = getGenomeAttribute('reduced_penetrance')
params.readcount_intervals             = getGenomeAttribute('readcount_intervals')
params.rtg_truthvcfs                   = getGenomeAttribute('rtg_truthvcfs')
params.sample_id_map                   = getGenomeAttribute('sample_id_map')
params.sequence_dictionary             = getGenomeAttribute('sequence_dictionary')
params.score_config_mt                 = getGenomeAttribute('score_config_mt')
params.score_config_snv                = getGenomeAttribute('score_config_snv')
params.score_config_sv                 = getGenomeAttribute('score_config_sv')
params.sdf                             = getGenomeAttribute('sdf')
params.svdb_query_bedpedbs             = getGenomeAttribute('svdb_query_bedpedbs')
params.svdb_query_dbs                  = getGenomeAttribute('svdb_query_dbs')
params.target_bed                      = getGenomeAttribute('target_bed')
params.variant_catalog                 = getGenomeAttribute('variant_catalog')
params.variant_consequences_snv        = getGenomeAttribute('variant_consequences_snv')
params.variant_consequences_sv         = getGenomeAttribute('variant_consequences_sv')
params.vep_filters                     = getGenomeAttribute('vep_filters')
params.vep_filters_scout_fmt           = getGenomeAttribute('vep_filters_scout_fmt')
params.vcf2cytosure_blacklist          = getGenomeAttribute('vcf2cytosure_blacklist')
params.vcfanno_resources               = getGenomeAttribute('vcfanno_resources')
params.vcfanno_toml                    = getGenomeAttribute('vcfanno_toml')
params.vcfanno_lua                     = getGenomeAttribute('vcfanno_lua')
params.vep_cache                       = getGenomeAttribute('vep_cache')
params.vep_plugin_files                = getGenomeAttribute('vep_plugin_files')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RAREDISEASE             } from './workflows/raredisease'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_raredisease_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_raredisease_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_RAREDISEASE {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    RAREDISEASE (
        samplesheet
    )

    emit:
    multiqc_report = RAREDISEASE.out.multiqc_report // channel: /path/to/multiqc_report.html

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_RAREDISEASE (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_RAREDISEASE.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//

def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}
