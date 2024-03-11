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
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RAREDISEASE             } from './workflows/raredisease'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_raredisease_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_raredisease_pipeline'

include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_raredisease_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

getGenomeAttribute('fasta','fasta')
getGenomeAttribute('fai','fai')
getGenomeAttribute('bwa','bwa')
getGenomeAttribute('bwamem2','bwamem2')
getGenomeAttribute('call_interval','call_interval')
getGenomeAttribute('cadd_resources','cadd_resources')
getGenomeAttribute('gcnvcaller_model','gcnvcaller_model')
getGenomeAttribute('gens_interval_list','gens_interval_list')
getGenomeAttribute('gens_pon_female','gens_pon_female')
getGenomeAttribute('gens_pon_male','gens_pon_male')
getGenomeAttribute('gens_gnomad_pos','gens_gnomad_pos')
getGenomeAttribute('gnomad_af','gnomad_af')
getGenomeAttribute('gnomad_af_idx','gnomad_af_idx')
getGenomeAttribute('intervals_wgs','intervals_wgs')
getGenomeAttribute('intervals_y','intervals_y')
getGenomeAttribute('known_dbsnp','known_dbsnp')
getGenomeAttribute('known_dbsnp_tbi','known_dbsnp_tbi')
getGenomeAttribute('mobile_element_references','mobile_element_references')
getGenomeAttribute('ml_model','ml_model')
getGenomeAttribute('mt_fasta','mt_fasta')
getGenomeAttribute('ploidy_model','ploidy_model')
getGenomeAttribute('reduced_penetrance','reduced_penetrance')
getGenomeAttribute('readcount_intervals','readcount_intervals')
getGenomeAttribute('rtg_truthvcfs','rtg_truthvcfs')
getGenomeAttribute('sequence_dictionary','sequence_dictionary')
getGenomeAttribute('score_config_mt','score_config_mt')
getGenomeAttribute('score_config_snv','score_config_snv')
getGenomeAttribute('score_config_sv','score_config_sv')
getGenomeAttribute('sdf','sdf')
getGenomeAttribute('svdb_query_dbs','svdb_query_dbs')
getGenomeAttribute('target_bed','target_bed')
getGenomeAttribute('variant_catalog','variant_catalog')
getGenomeAttribute('variant_consequences_snv','variant_consequences_snv')
getGenomeAttribute('variant_consequences_sv','variant_consequences_sv')
getGenomeAttribute('vep_filters','vep_filters')
getGenomeAttribute('vep_filters_scout_fmt','vep_filters_scout_fmt')
getGenomeAttribute('vcf2cytosure_blacklist','vcf2cytosure_blacklist')
getGenomeAttribute('vcfanno_resources','vcfanno_resources')
getGenomeAttribute('vcfanno_toml','vcfanno_toml')
getGenomeAttribute('vcfanno_lua','vcfanno_lua')
getGenomeAttribute('vep_cache','vep_cache')
getGenomeAttribute('vep_plugin_files','vep_plugin_files')
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
