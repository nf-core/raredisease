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
    fastq
    alignment
    samples
    case_info

    main:

    //
    // WORKFLOW: Run pipeline
    //

    if (params.skip_tools) {
        skip_eklipse             = params.skip_tools.split(',').contains('eklipse')
        skip_fastp               = params.skip_tools.split(',').contains('fastp')
        skip_fastqc              = params.skip_tools.split(',').contains('fastqc')
        skip_gens                = params.skip_tools.split(',').contains('gens')
        skip_germlinecnvcaller   = params.skip_tools.split(',').contains('germlinecnvcaller')
        skip_haplogrep3          = params.skip_tools.split(',').contains('haplogrep3')
        skip_ngsbits             = params.skip_tools.split(',').contains('ngsbits')
        skip_peddy               = params.skip_tools.split(',').contains('peddy')
        skip_qualimap            = params.skip_tools.split(',').contains('qualimap')
        skip_smncopynumbercaller = params.skip_tools.split(',').contains('smncopynumbercaller')
        skip_vcf2cytosure        = params.skip_tools.split(',').contains('vcf2cytosure')
    }

    if (params.skip_subworkflows) {
        skip_me_annotation         = params.skip_tools.split(',').contains('me_annotation')
        skip_me_calling            = params.skip_tools.split(',').contains('me_calling')
        skip_mt_annotation         = params.skip_tools.split(',').contains('mt_annotation')
        skip_mt_subsample          = params.skip_tools.split(',').contains('mt_subsample')
        skip_repeat_annotation     = params.skip_tools.split(',').contains('repeat_annotation')
        skip_repeat_calling        = params.skip_tools.split(',').contains('repeat_calling')
        skip_snv_annotation        = params.skip_tools.split(',').contains('snv_annotation')
        skip_snv_calling           = params.skip_tools.split(',').contains('snv_calling')
        skip_sv_annotation         = params.skip_tools.split(',').contains('sv_annotation')
        skip_sv_calling            = params.skip_tools.split(',').contains('sv_calling')
        skip_generate_clinical_set = params.skip_tools.split(',').contains('generate_clinical_set')
    }

    RAREDISEASE (
        fastq,
        alignment,
        samples,
        case_info,
        params.analysis_type,
        skip_me_calling,
        skip_me_annotation,
        skip_mt_annotation,
        skip_mt_subsample,
        skip_repeat_annotation,
        skip_repeat_calling,
        skip_snv_annotation,
        skip_snv_calling,
        skip_sv_annotation,
        skip_sv_calling,
        skip_generate_clinical_set,
        skip_eklipse,
        skip_fastp,
        skip_fastqc,
        skip_gens,
        skip_germlinecnvcaller,
        skip_haplogrep3,
        skip_ngsbits,
        skip_peddy,
        skip_qualimap,
        skip_smncopynumbercaller,
        skip_vcf2cytosure
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
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )
    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_RAREDISEASE (
        PIPELINE_INITIALISATION.out.reads,
        PIPELINE_INITIALISATION.out.align,
        PIPELINE_INITIALISATION.out.samples,
        PIPELINE_INITIALISATION.out.case_info
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
