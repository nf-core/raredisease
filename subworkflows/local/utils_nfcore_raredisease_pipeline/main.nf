//
// Subworkflow with functionality specific to the nf-core/raredisease pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .tap { ch_original_input }
        .map { meta, fastq1, fastq2, spring1, spring2, bam, bai -> meta.id }
        .reduce([:]) { counts, sample -> //get counts of each sample in the samplesheet - for groupTuple
            counts[sample] = (counts[sample] ?: 0) + 1
            counts
        }
        .combine( ch_original_input )
        .map { counts, meta, fastq1, fastq2, spring1, spring2, bam, bai ->
            def new_meta = meta + [num_lanes:counts[meta.id]]
            if (fastq1 && fastq2) {
                new_meta += [read_group: generateReadGroupLine(fastq1, meta, params)]
                return [new_meta + [single_end: false, data_type: "fastq_gz"], [fastq1, fastq2]]
            } else if (fastq1 && !fastq2) {
                new_meta += [read_group: generateReadGroupLine(fastq1, meta, params)]
                return [new_meta + [single_end: true, data_type: "fastq_gz"], [fastq1]]
            } else if (spring1 && spring2) {
                new_meta += [read_group: generateReadGroupLine(spring1, meta, params)]
                return [new_meta + [single_end: false, data_type: "separate_spring"], [spring1, spring2]]
            } else if (spring1 && !spring2) {
                new_meta += [read_group: generateReadGroupLine(spring1, meta, params)]
                return [new_meta + [single_end: false, data_type: "interleaved_spring"], [spring1]]
            } else if (bam && bai) {
                new_meta += [read_group: generateReadGroupLine(bam, meta, params)]
                return [new_meta, [bam, bai]]
            }
        }
        .tap{ ch_input_counts }
        .map { meta, files -> files }
        .reduce([:]) { counts, files -> //get line number for each row to construct unique sample ids
            counts[files] = counts.size() + 1
            return counts
        }
        .combine( ch_input_counts )
        .map { lineno, meta, files -> //append line number to sampleid
            def new_meta = meta + [id:meta.id+"_LNUMBER"+lineno[files]]
            return [ new_meta, files ]
        }
        .tap { ch_samplesheet }
        .branch { meta, files  ->
            fastq: !files[0].toString().endsWith("bam")
                return [meta, files]
            align: files[0].toString().endsWith("bam")
                return [meta, files]
        }
        .set {ch_samplesheet_by_type}

    ch_samples  = ch_samplesheet.map { meta, files ->
                    def new_id = meta.sample
                    def new_meta = meta - meta.subMap('lane', 'read_group') + [id:new_id]
                    return new_meta
                    }.unique()

    ch_case_info = ch_samples.toList().map { createCaseChannel(it) }

    emit:
    reads     = ch_samplesheet_by_type.fastq
    align     = ch_samplesheet_by_type.align
    samples   = ch_samples
    case_info = ch_case_info
    versions  = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def generateReadGroupLine(file, meta, params) {
    return "\'@RG\\tID:" + file.simpleName + "_" + meta.lane + "\\tPL:" + params.platform.toUpperCase() + "\\tSM:" + meta.id + "\'"
}

def boolean isNonZeroNonEmpty(value) {
        return (value instanceof String && value != "" && value != "0") ||
                (value instanceof Number && value != 0)
    }

    // Function to get a list of metadata (e.g. case id) for the case [ meta ]
def createCaseChannel(List rows) {
    def case_info    = [:]
    def probands     = [] as Set
    def upd_children = [] as Set
    def father       = ""
    def mother       = ""

    rows.each { item ->
        if (item?.phenotype == 2) {
            probands << item.sample
        }
        if (isNonZeroNonEmpty(item?.paternal) && isNonZeroNonEmpty(item?.maternal)) {
            upd_children << item.sample
        }
        if (isNonZeroNonEmpty(item?.paternal)) {
            father = item.paternal
        }
        if (isNonZeroNonEmpty(item?.maternal)) {
            mother = item.maternal
        }
    }

    case_info.father       = father
    case_info.mother       = mother
    case_info.probands     = probands.toList()
    case_info.upd_children = upd_children.toList()
    case_info.id           = rows[0].case_id

    return case_info
}

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate parameters
//

def checkRequiredParameters(params) {
    def mandatoryParams = [
        "analysis_type",
        "fasta",
        "input",
        "intervals_wgs",
        "intervals_y",
        "variant_caller"
    ]

    // Static requirements that are not influenced by user-defined skips
    def staticRequirements   = [
        analysis_type_wes        : ["target_bed"],
        variant_caller_sentieon  : ["ml_model"],
        run_rtgvcfeval           : ["rtg_truthvcfs"]
    ]

    // Requirements that can be modified by the user using either skip_tools or skip_subworkflows here
    def dynamicRequirements = [
        repeat_calling           : ["variant_catalog"],
        repeat_annotation        : ["variant_catalog"],
        snv_calling              : ["genome"],
        snv_annotation           : ["genome", "vcfanno_resources", "vcfanno_toml", "vep_cache", "vep_cache_version",
                                    "gnomad_af", "score_config_snv", "variant_consequences_snv"],
        sv_annotation            : ["genome", "vep_cache", "vep_cache_version", "score_config_sv", "variant_consequences_sv"],
        mt_annotation            : ["genome", "mito_name", "vcfanno_resources", "vcfanno_toml", "score_config_mt",
                                    "vep_cache_version", "vep_cache", "variant_consequences_snv"],
        me_calling               : ["mobile_element_references"],
        me_annotation            : ["mobile_element_svdb_annotations", "variant_consequences_snv"],
        gens                     : ["gens_gnomad_pos", "gens_interval_list", "gens_pon_female", "gens_pon_male"],
        germlinecnvcaller        : ["ploidy_model", "gcnvcaller_model", "readcount_intervals"],
        smncopynumbercaller      : ["genome"]
    ]

    def missingParamsCount = 0

    staticRequirements.each { condition, paramsList ->
        if ((condition == "analysis_type_wes" && params.analysis_type == "wes") ||
            (condition == "variant_caller_sentieon" && params.variant_caller == "sentieon") ||
            (condition == "run_rtgvcfeval" && params.run_rtgvcfeval)) {
                mandatoryParams += paramsList
        }
    }

    all_skips = params.skip_subworkflows+","+params.skip_tools
    dynamicRequirements.each { condition, paramsList ->
        if (!all_skips.split(',').contains(condition)) {
                mandatoryParams += paramsList
        }
    }

    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('sv_annotation')) && !params.svdb_query_bedpedbs && !params.svdb_query_dbs) {
        println("params.svdb_query_bedpedbs or params.svdb_query_dbs should be set.")
        missingParamsCount += 1
    }

    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('generate_clinical_set')) ) {
        if (!params.vep_filters && !params.vep_filters_scout_fmt) {
            println("params.vep_filters or params.vep_filters_scout_fmt should be set.")
            missingParamsCount += 1
        } else if (params.vep_filters && params.vep_filters_scout_fmt) {
            println("Either params.vep_filters or params.vep_filters_scout_fmt should be set.")
            missingParamsCount += 1
        }
    }

    mandatoryParams.unique().each { param ->
        if (params[param] == null) {
            println("params." + param + " not set.")
            missingParamsCount += 1
        }
    }

    if (missingParamsCount > 0) {
        error("\nSet missing parameters and restart the run. For more information please check usage documentation on github.")
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

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

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {

    def align_text           = []

    align_text = [
        params.aligner.equals("bwa")      ? "BWA (Li, 2013),"                        :"",
        params.aligner.equals("bwamem2")  ? "BWA-MEM2 (Vasimuddin et al., 2019),"    : "",
        params.aligner.equals("bwameme")  ? "BWA-MEME (Jung et al., 2022),"          : "",
        params.aligner.equals("sentieon") ? "Sentieon DNASeq (Kendig et al., 2019)," : "",
        params.aligner.equals("sentieon") ? "Sentieon Tools (Freed et al., 2017),"   : ""
    ]

    def concat_text = align_text

    def citation_text = [ "Tools used in the workflow included:" ] + concat_text.unique(false) { a, b -> a <=> b } - ""
    return citation_text.join(' ').trim()
}

def toolBibliographyText() {

    def align_text           = []

    align_text = [
        params.aligner.equals("bwa") ? "<li>Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM (arXiv:1303.3997). arXiv. http://arxiv.org/abs/1303.3997</li>" :"",
        params.aligner.equals("bwamem2") ? "<li>Vasimuddin, Md., Misra, S., Li, H., & Aluru, S. (2019). Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. 2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS), 314–324. https://doi.org/10.1109/IPDPS.2019.00041</li>" : "",
        params.aligner.equals("bwameme") ? "<li>Jung Y, Han D. BWA-MEME: BWA-MEM emulated with a machine learning approach. Bioinformatics. 2022;38(9):2404-2413. doi:10.1093/bioinformatics/btac137</li>" : "",
        params.aligner.equals("sentieon") ? "<li>Kendig, K. I., Baheti, S., Bockol, M. A., Drucker, T. M., Hart, S. N., Heldenbrand, J. R., Hernaez, M., Hudson, M. E., Kalmbach, M. T., Klee, E. W., Mattson, N. R., Ross, C. A., Taschuk, M., Wieben, E. D., Wiepert, M., Wildman, D. E., & Mainzer, L. S. (2019). Sentieon DNASeq Variant Calling Workflow Demonstrates Strong Computational Performance and Accuracy. Frontiers in Genetics, 10, 736. https://doi.org/10.3389/fgene.2019.00736</li>" : "",
        params.aligner.equals("sentieon") ? "<li>Freed, D., Aldana, R., Weber, J. A., & Edwards, J. S. (2017). The Sentieon Genomics Tools—A fast and accurate solution to variant calling from next-generation sequence data (p. 115717). bioRxiv. https://doi.org/10.1101/115717</li>" : ""
    ]

    def concat_text = align_text 
    def reference_text = concat_text.unique(false) { a, b -> a <=> b } - ""
    return reference_text.join(' ').trim()
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
