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
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
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
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    monochrome_logs   // boolean: Disable ANSI colour codes in log output
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

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

    def before_text = ""
    def after_text = ""
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/raredisease ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/raredisease/blob/master/CITATIONS.md
"""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
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
    checkRequiredParameters(params)

    //
    // Create channel from input file provided through params.input
    //
    channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .tap { ch_original_input }
        .map { meta, _fastq1, _fastq2, _spring1, _spring2, _bam, _bai, _cram, _crai, _vcf, _tbi, _type -> meta.id }
        .reduce([:]) { counts, sample -> //get counts of each sample in the samplesheet - for groupTuple
            counts[sample] = (counts[sample] ?: 0) + 1
            counts
        }
        .combine( ch_original_input )
        .map { counts, meta, fastq1, fastq2, spring1, spring2, bam, bai, cram, crai, vcf, tbi, type ->
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
                return [new_meta + [data_type: "bam"], [bam, bai]]
            } else if (cram && crai) {
                new_meta += [read_group: generateReadGroupLine(cram, meta, params)]
                return [new_meta + [data_type: "cram"], [cram, crai]]
            } else if (vcf && tbi && type) {
                return [new_meta + [data_type: "${type}_vcf"], [vcf, tbi]]
            }
        }
        .tap{ ch_input_counts }
        .map { _meta, files -> files }
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
            fastq:     meta.data_type in ["fastq_gz", "separate_spring", "interleaved_spring"]
                return [meta, files]
            align:     meta.data_type in ["bam", "cram"]
                return [meta, files]
            precalled: meta.data_type.endsWith("_vcf")
                return [meta, files]
        }
        .set {ch_samplesheet_by_type}

    ch_samples  = ch_samplesheet.map { meta, _files ->
                    def new_id = meta.sample
                    def new_meta = meta - meta.subMap('lane', 'read_group') + [id:new_id]
                    return new_meta
                    }.unique()

    ch_case_info = ch_samples.toList().map { it -> validateNoMixedCaseInput(it); createCaseChannel(it) }

    ch_precalled_vcfs = ch_samplesheet_by_type.precalled
        .toList()
        .map { rows -> extractPrecalledVcfs(rows) }

    emit:
    reads          = ch_samplesheet_by_type.fastq // channel: [ val(meta), [ path(reads) ] ]
    align          = ch_samplesheet_by_type.align // channel: [ val(meta), [ path(bam/cram), path(bai/crai) ] ]
    samples        = ch_samples                   // channel: [ val(meta) ]
    case_info      = ch_case_info                 // channel: [ val(case_info) ]
    precalled_vcfs = ch_precalled_vcfs             // channel: [ val([snv:[vcf,tbi]|null, sv:[...]|null, mt:[...]|null]) ]
    versions       = ch_versions                  // channel: [ path(versions) ]
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

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/**
 * Creates a channel from a file path if provided, otherwise returns a fallback channel
 * @param filePath The path to the file (can be null)
 * @param valueFallback If true, returns channel.value([]) when filePath is null; otherwise returns channel.empty() (default: false)
 * @return Channel with collected file path or fallback channel
 */
def channelFromPath(filePath, valueFallback = false) {
    if (!filePath) {
        return valueFallback ? channel.value([]) : channel.empty()
    }
    return channel.fromPath(filePath).collect()
}

/**
 * Creates a channel from a file path, maps it to [id, file] format, and collects
 * @param filePath The path to the file (can be null)
 * @param doubleEmpty If true, returns channel.value([[:], []]) when filePath is null; otherwise returns channel.empty() (default: false)
 * @param customId The custom ID to be used in meta.id (default: null)
 * @return Channel with [[id:name], file] format and collected, or fallback channel
 */
def channelFromPathWithMeta(filePath, doubleEmpty = false, customId = null) {
    if (!filePath) {
        return doubleEmpty ? channel.value([[:], []]) : channel.empty()
    }
    return channel.fromPath(filePath).map { file ->
        def meta_id = customId ?: file.simpleName
        return [[id: meta_id], file]
    }.collect()
}

/**
 * Creates a channel from a samplesheet file using samplesheetToList, or returns a fallback channel
 * @param samplesheetPath The path to the samplesheet file (can be null)
 * @param schemaPath The path to the JSON schema file for validation
 * @param collect If true, calls .collect() on the channel (default: true)
 * @return Channel from samplesheet list or channel.empty()
 */
def channelFromSamplesheet(samplesheetPath, schemaPath, collect = true) {
    if (!samplesheetPath) {
        return channel.empty()
    }
    def ch_out = channel.fromList(samplesheetToList(samplesheetPath, schemaPath))
    return collect ? ch_out.collect() : ch_out
}

def boolean hasSpringInput() {
    return file(params.input).readLines().any { line -> line.contains('.spring') }
}

// Checks whether any samplesheet row has its 'type' column set to the given precalled-VCF type (snv/sv/mt)
def boolean hasPrecalledVcfOfType(String type) {
    def lines = file(params.input).readLines()
    if (!lines) {
        return false
    }
    def header   = lines[0].split(',', -1)*.trim()
    def type_idx = header.indexOf('type')
    if (type_idx == -1) {
        return false
    }
    return lines.drop(1).any { line ->
        def fields = line.split(',', -1)
        type_idx < fields.size() && fields[type_idx].trim() == type
    }
}

// Per-type wrappers used to compute the has_precalled_* emits and gate calling/annotation downstream
def boolean hasPrecalledSnvVcf() {
    return hasPrecalledVcfOfType('snv')
}

def boolean hasPrecalledSvVcf() {
    return hasPrecalledVcfOfType('sv')
}

def boolean hasPrecalledMtVcf() {
    return hasPrecalledVcfOfType('mt')
}

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

// A case must be either fully precalled (vcf/tbi/type rows) or fully processed from raw/aligned reads, never both
def validateNoMixedCaseInput(List rows) {
    def groups_by_case = [:]
    rows.each { row ->
        def group = row.data_type.endsWith('_vcf') ? 'vcf' : 'align'
        groups_by_case.computeIfAbsent(row.case_id) { [] as Set }.add(group)
    }
    groups_by_case.each { case_id, groups ->
        if (groups.size() > 1) {
            error("Case '${case_id}' mixes precalled VCF input (vcf/tbi/type columns) with fastq/spring/bam/cram input. A case must be either fully precalled or fully processed from raw/aligned reads, not both.")
        }
    }
}

// Function to collect precalled vcf/tbi pairs per variant type from rows tagged with a "*_vcf" data_type
def extractPrecalledVcfs(List rows) {
    def precalled = [snv: null, sv: null, mt: null]
    rows.each { meta, files ->
        def type = meta.data_type - "_vcf"
        if (precalled[type] && precalled[type] != files) {
            error("Conflicting precalled '${type}' VCFs supplied in samplesheet for the same case.")
        }
        precalled[type] = files
    }
    return precalled
}

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
    validatePrecalledVcfCoverage()
}

// A case with any precalled VCF has no fastq/bam/cram rows (enforced by validateNoMixedCaseInput), so every
// relevant type must either have a precalled VCF or have its calling explicitly skipped, otherwise it would
// silently run calling with no input data
def validatePrecalledVcfCoverage() {
    def has_snv = hasPrecalledSnvVcf()
    def has_sv  = hasPrecalledSvVcf()
    def has_mt  = hasPrecalledMtVcf()
    if (!has_snv && !has_sv && !has_mt) {
        return
    }
    def run_mt  = params.analysis_type.matches("wgs|mito") || params.run_mt_for_wes
    def missing = []
    if (!has_snv && !parseSkipList(params.skip_subworkflows, 'snv_calling')) {
        missing << 'snv'
    }
    if (!has_sv && !parseSkipList(params.skip_subworkflows, 'sv_calling')) {
        missing << 'sv'
    }
    if (run_mt && !has_mt && !parseSkipList(params.skip_subworkflows, 'mt_calling')) {
        missing << 'mt'
    }
    if (missing) {
        error("The samplesheet supplies a precalled VCF for at least one variant type, making this a fully-precalled case with no fastq/bam/cram input for calling. But no precalled VCF was supplied for: ${missing.join(', ')}. Either add a precalled VCF for ${missing.join(', ')}, or skip that calling explicitly via --skip_subworkflows.")
    }
}


//
// Initialize skip parameters
//
def parseSkipList(paramValue, toolName) {
    return paramValue ? paramValue.split(',').contains(toolName) : false
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

    def all_skips = params.skip_subworkflows+","+params.skip_tools
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
        } else {
            def filtersFile = params.vep_filters_scout_fmt ?: params.vep_filters
            def nonHeaderLines = file(filtersFile).readLines().count { line -> !line.startsWith('#') && line.trim() }
            if (nonHeaderLines == 0) {
                error("The file '${filtersFile}' contains no records (only headers or empty lines). The clinical set will contain 0 variants.")
            }
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
    def repeats_text         = []
    def snv_annotation_text  = []
    def snv_calls_text       = []
    def sv_annotation_text   = []
    def sv_calls_text        = []
    def mt_annotation_text   = []
    def qc_bam_text          = []
    def me_calls_text        = []
    def me_annotation_text   = []
    def preprocessing_text   = []
    def other_citation_text  = []

    align_text = [
        params.aligner.equals("bwa")      ? "BWA (Li, 2013),"                        :"",
        params.aligner.equals("bwamem2")  ? "BWA-MEM2 (Vasimuddin et al., 2019),"    : "",
        params.aligner.equals("bwameme")  ? "BWA-MEME (Jung et al., 2022),"          : "",
        params.aligner.equals("sentieon") ? "Sentieon DNASeq (Kendig et al., 2019)," : "",
        params.aligner.equals("sentieon") ? "Sentieon Tools (Freed et al., 2017),"   : ""
    ]
    repeats_text = [
        (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('repeat_calling')) && params.analysis_type.equals("wgs"))    ? "ExpansionHunter (Dolzhenko et al., 2019)," : "",
        (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('repeat_annotation')) && params.analysis_type.equals("wgs")) ? "stranger (Nilsson & Magnusson, 2021)," : ""
    ]
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('snv_annotation'))) {
        snv_annotation_text = [
            "CADD (Rentzsch et al., 2019, 2021),",
            "Vcfanno (Pedersen et al., 2016),",
            "VEP (McLaren et al., 2016),",
            "Genmod (Magnusson et al., 2018),"
        ]
    }
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('snv_calling'))) {
        snv_calls_text = [
            params.variant_caller.equals("deepvariant") ? "DeepVariant (Poplin et al., 2018),"      : "",
            params.variant_caller.equals("sentieon")    ? "Sentieon DNAscope (Freed et al., 2022)," : "",
            params.run_mt_for_wes ? "Haplocheck (Weissensteiner et al., 2021)," : "",
            "GLnexus (Yun et al., 2021),"
        ]
    }
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('sv_annotation'))) {
        sv_annotation_text = [
            "SVDB (Eisfeldt et al., 2017),",
            "VEP (McLaren et al., 2016),",
            "Genmod (Magnusson et al., 2018),"
        ]
    }
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('sv_calling'))) {
        sv_calls_text = [
            params.analysis_type.equals("wgs") ? "CNVnator (Abyzov et al., 2011)," : "",
            params.analysis_type.equals("wgs") ? "TIDDIT (Eisfeldt et al., 2017)," : "",
            "Manta (Chen et al., 2016),"
        ]
    }
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('mt_annotation')) && (params.analysis_type.equals("wgs") || params.run_mt_for_wes)) {
        mt_annotation_text = [
            "CADD (Rentzsch et al., 2019, 2021),",
            "VEP (McLaren et al., 2016),",
            "Vcfanno (Pedersen et al., 2016),",
            "Genmod (Magnusson et al., 2018),"
        ]
    }
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('me_annotation')) && params.analysis_type.equals("wgs")) {
        me_annotation_text = [
            "VEP (McLaren et al., 2016),",
            "SVDB (Eisfeldt et al., 2017),"
        ]
    }
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('me_calling')) && params.analysis_type.equals("wgs")) {
        me_calls_text = [
            "SVDB (Eisfeldt et al., 2017),",
            "RetroSeq (Keane et al., 2013),"
        ]
    }
    qc_bam_text = [
        "Picard (Broad Institute, 2023)",
        "Sambamba (Tarasov et al., 2015),",
        "TIDDIT (Eisfeldt et al., 2017),",
        "UCSC Bigwig and Bigbed (Kent et al., 2010),",
        (params.verifybamid_svd_bed && params.verifybamid_svd_mu && params.verifybamid_svd_ud) ? "VerifyBamID2 (Zhang et al., 2020)," : "",
        "Mosdepth (Pedersen & Quinlan, 2018),"
    ]
    preprocessing_text = [
        "FastQC (Andrews 2010),",
        (params.skip_tools && params.skip_tools.split(',').contains('fastp')) ? "" : "Fastp (Chen, 2023),",
        hasSpringInput() ? "Spring (Chandak et al., 2019)," : ""
    ]
    other_citation_text = [
        "BCFtools (Danecek et al., 2021),",
        "BEDTools (Quinlan & Hall, 2010),",
        "GATK (McKenna et al., 2010),",
        "MultiQC (Ewels et al. 2016),",
        (params.skip_tools && params.skip_tools.split(',').contains('peddy')) ? "" : "Peddy (Pedersen & Quinlan, 2017),",
        params.run_rtgvcfeval ? "RTG Tools (Cleary et al., 2015)," : "",
        "SAMtools (Li et al., 2009),",
        (!(params.skip_tools && params.skip_tools.split(',').contains('smncopynumbercaller')) && params.analysis_type.equals("wgs")) ? "SMNCopyNumberCaller (Chen et al., 2020)," : "",
        "Tabix (Li, 2011)",
        "."
    ]

    def concat_text = align_text +
                        repeats_text         +
                        snv_annotation_text  +
                        snv_calls_text       +
                        sv_annotation_text   +
                        sv_calls_text        +
                        mt_annotation_text   +
                        qc_bam_text          +
                        me_calls_text        +
                        me_annotation_text   +
                        preprocessing_text   +
                        other_citation_text

    def citation_text = [ "Tools used in the workflow included:" ] + concat_text.unique(false) { a, b -> a <=> b } - ""
    return citation_text.join(' ').trim()
}

def toolBibliographyText() {

    def align_text           = []
    def repeats_text         = []
    def snv_annotation_text  = []
    def snv_calls_text       = []
    def sv_annotation_text   = []
    def sv_calls_text        = []
    def mt_annotation_text   = []
    def qc_bam_text          = []
    def me_calls_text        = []
    def me_annotation_text   = []
    def preprocessing_text   = []
    def other_citation_text  = []

    align_text = [
        params.aligner.equals("bwa") || params.mt_aligner.equals("bwa") ? "<li>Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM (arXiv:1303.3997). arXiv. http://arxiv.org/abs/1303.3997</li>" :"",
        params.aligner.equals("bwamem2") || params.mt_aligner.equals("bwamem2") ? "<li>Vasimuddin, Md., Misra, S., Li, H., & Aluru, S. (2019). Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. 2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS), 314–324. https://doi.org/10.1109/IPDPS.2019.00041</li>" : "",
        params.aligner.equals("bwameme") ? "<li>Jung Y, Han D. BWA-MEME: BWA-MEM emulated with a machine learning approach. Bioinformatics. 2022;38(9):2404-2413. doi:10.1093/bioinformatics/btac137</li>" : "",
        params.aligner.equals("sentieon") || params.mt_aligner.equals("sentieon") ? "<li>Kendig, K. I., Baheti, S., Bockol, M. A., Drucker, T. M., Hart, S. N., Heldenbrand, J. R., Hernaez, M., Hudson, M. E., Kalmbach, M. T., Klee, E. W., Mattson, N. R., Ross, C. A., Taschuk, M., Wieben, E. D., Wiepert, M., Wildman, D. E., & Mainzer, L. S. (2019). Sentieon DNASeq Variant Calling Workflow Demonstrates Strong Computational Performance and Accuracy. Frontiers in Genetics, 10, 736. https://doi.org/10.3389/fgene.2019.00736</li>" : "",
        params.aligner.equals("sentieon") || params.mt_aligner.equals("sentieon")? "<li>Freed, D., Aldana, R., Weber, J. A., & Edwards, J. S. (2017). The Sentieon Genomics Tools—A fast and accurate solution to variant calling from next-generation sequence data (p. 115717). bioRxiv. https://doi.org/10.1101/115717</li>" : ""
    ]
    repeats_text = [
        (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('repeat_calling')) && params.analysis_type.equals("wgs") )    ? "<li>Dolzhenko, E., Deshpande, V., Schlesinger, F., Krusche, P., Petrovski, R., Chen, S., Emig-Agius, D., Gross, A., Narzisi, G., Bowman, B., Scheffler, K., van Vugt, J. J. F. A., French, C., Sanchis-Juan, A., Ibáñez, K., Tucci, A., Lajoie, B. R., Veldink, J. H., Raymond, F. L., … Eberle, M. A. (2019). ExpansionHunter: A sequence-graph-based tool to analyze variation in short tandem repeat regions. Bioinformatics, 35(22), 4754–4756. https://doi.org/10.1093/bioinformatics/btz431</li>" : "",
        (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('repeat_annotation')) && params.analysis_type.equals("wgs") ) ? "<li>Nilsson, D., & Magnusson, M. (2021). Moonso/stranger v0.7.1 (v0.7.1) [Computer software]. Zenodo. https://doi.org/10.5281/ZENODO.4548873</li>" : ""
    ]
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('snv_annotation'))) {
        snv_annotation_text = [
            "<li>Rentzsch, P., Schubach, M., Shendure, J., & Kircher, M. (2021). CADD-Splice—Improving genome-wide variant effect prediction using deep learning-derived splice scores. Genome Medicine, 13(1), 31. https://doi.org/10.1186/s13073-021-00835-9</li>",
            "<li>Rentzsch, P., Witten, D., Cooper, G. M., Shendure, J., & Kircher, M. (2019). CADD: Predicting the deleteriousness of variants throughout the human genome. Nucleic Acids Research, 47(D1), D886–D894. https://doi.org/10.1093/nar/gky1016</li>",
            "<li>Pedersen, B. S., Layer, R. M., & Quinlan, A. R. (2016). Vcfanno: Fast, flexible annotation of genetic variants. Genome Biology, 17(1), 118. https://doi.org/10.1186/s13059-016-0973-5</li>",
            "<li>McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R. S., Thormann, A., Flicek, P., & Cunningham, F. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1), 122. https://doi.org/10.1186/s13059-016-0974-4</li>",
            "<li>Magnusson, M., Hughes, T., Glabilloy, & Bitdeli Chef. (2018). genmod: Version 3.7.3 (3.7.3) [Computer software]. Zenodo. https://doi.org/10.5281/ZENODO.3841142</li>"
        ]
    }
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('snv_calling'))) {
        snv_calls_text = [
            params.variant_caller.equals("deepvariant") ? "<li>Poplin, R., Chang, P.-C., Alexander, D., Schwartz, S., Colthurst, T., Ku, A., Newburger, D., Dijamco, J., Nguyen, N., Afshar, P. T., Gross, S. S., Dorfman, L., McLean, C. Y., & DePristo, M. A. (2018). A universal SNP and small-indel variant caller using deep neural networks. Nature Biotechnology, 36(10), 983–987. https://doi.org/10.1038/nbt.4235</li>" : "",
            params.variant_caller.equals("sentieon") ? "<li>Freed, D., Pan, R., Chen, H., Li, Z., Hu, J., & Aldana, R. (2022). DNAscope: High accuracy small variant calling using machine learning [Preprint]. Bioinformatics. https://doi.org/10.1101/2022.05.20.492556</li>" : "",
            params.run_mt_for_wes ? "<li>Weissensteiner, H., Forer, L., Fendt, L., Kheirkhah, A., Salas, A., Kronenberg, F., & Schoenherr, S. (2021). Contamination detection in sequencing studies using the mitochondrial phylogeny. Genome Research, 31(2), 309–316. https://doi.org/10.1101/gr.256545.119</li>" : "",
            "<li>Yun, T., Li, H., Chang, P.-C., Lin, M. F., Carroll, A., & McLean, C. Y. (2021). Accurate, scalable cohort variant calls using DeepVariant and GLnexus. Bioinformatics, 36(24), 5582–5589. https://doi.org/10.1093/bioinformatics/btaa1081</li>"
        ]
    }

    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('sv_annotation'))) {
        sv_annotation_text = [
            "<li>Eisfeldt, J., Vezzi, F., Olason, P., Nilsson, D., & Lindstrand, A. (2017). TIDDIT, an efficient and comprehensive structural variant caller for massive parallel sequencing data. F1000Research, 6, 664. https://doi.org/10.12688/f1000research.11168.2</li>",
            "<li>McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R. S., Thormann, A., Flicek, P., & Cunningham, F. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1), 122. https://doi.org/10.1186/s13059-016-0974-4</li>",
            "<li>Magnusson, M., Hughes, T., Glabilloy, & Bitdeli Chef. (2018). genmod: Version 3.7.3 (3.7.3) [Computer software]. Zenodo. https://doi.org/10.5281/ZENODO.3841142</li>"
        ]
    }
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('sv_calling'))) {
        sv_calls_text = [
            params.analysis_type.equals("wgs") ? "<li>Abyzov, A., Urban, A. E., Snyder, M., & Gerstein, M. (2011). CNVnator: An approach to discover, genotype, and characterize typical and atypical CNVs from family and population genome sequencing. Genome Research, 21(6), 974–984. https://doi.org/10.1101/gr.114876.110</li>" : "",
            params.analysis_type.equals("wgs") ? "<li>Eisfeldt, J., Vezzi, F., Olason, P., Nilsson, D., & Lindstrand, A. (2017). TIDDIT, an efficient and comprehensive structural variant caller for massive parallel sequencing data. F1000Research, 6, 664. https://doi.org/10.12688/f1000research.11168.2</li>" : "",
            "<li>Chen, X., Schulz-Trieglaff, O., Shaw, R., Barnes, B., Schlesinger, F., Källberg, M., Cox, A. J., Kruglyak, S., & Saunders, C. T. (2016). Manta: Rapid detection of structural variants and indels for germline and cancer sequencing applications. Bioinformatics, 32(8), 1220–1222. https://doi.org/10.1093/bioinformatics/btv710</li>"
        ]
    }

    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('mt_annotation')) && (params.analysis_type.equals("wgs") || params.run_mt_for_wes)) {
        mt_annotation_text = [
            "<li>Rentzsch, P., Schubach, M., Shendure, J., & Kircher, M. (2021). CADD-Splice—Improving genome-wide variant effect prediction using deep learning-derived splice scores. Genome Medicine, 13(1), 31. https://doi.org/10.1186/s13073-021-00835-9</li>",
            "<li>Rentzsch, P., Witten, D., Cooper, G. M., Shendure, J., & Kircher, M. (2019). CADD: Predicting the deleteriousness of variants throughout the human genome. Nucleic Acids Research, 47(D1), D886–D894. https://doi.org/10.1093/nar/gky1016</li>",
            "<li>Pedersen, B. S., Layer, R. M., & Quinlan, A. R. (2016). Vcfanno: Fast, flexible annotation of genetic variants. Genome Biology, 17(1), 118. https://doi.org/10.1186/s13059-016-0973-5</li>",
            "<li>McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R. S., Thormann, A., Flicek, P., & Cunningham, F. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1), 122. https://doi.org/10.1186/s13059-016-0974-4</li>",
            "<li>Magnusson, M., Hughes, T., Glabilloy, & Bitdeli Chef. (2018). genmod: Version 3.7.3 (3.7.3) [Computer software]. Zenodo. https://doi.org/10.5281/ZENODO.3841142</li>"
        ]
    }
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('me_annotation')) && params.analysis_type.equals("wgs")) {
        me_annotation_text = [
            "<li>McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R. S., Thormann, A., Flicek, P., & Cunningham, F. (2016). The Ensembl Variant Effect Predictor. Genome Biology, 17(1), 122. https://doi.org/10.1186/s13059-016-0974-4</li>",
            "<li>Eisfeldt, J., Vezzi, F., Olason, P., Nilsson, D., & Lindstrand, A. (2017). TIDDIT, an efficient and comprehensive structural variant caller for massive parallel sequencing data. F1000Research, 6, 664. https://doi.org/10.12688/f1000research.11168.2</li>"
        ]
    }
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('me_calling')) && params.analysis_type.equals("wgs")) {
        me_calls_text = [
            "<li>Eisfeldt, J., Vezzi, F., Olason, P., Nilsson, D., & Lindstrand, A. (2017). TIDDIT, an efficient and comprehensive structural variant caller for massive parallel sequencing data. F1000Research, 6, 664. https://doi.org/10.12688/f1000research.11168.2</li>",
            "<li>Keane, T. M., Wong, K., & Adams, D. J. (2013). RetroSeq: Transposable element discovery from next-generation sequencing data. Bioinformatics, 29(3), 389–390. https://doi.org/10.1093/bioinformatics/bts697</li>"
        ]
    }
    qc_bam_text = [
        "<li>Broad Institute. (2023). Picard Tools. In Broad Institute, GitHub repository. http://broadinstitute.github.io/picard/</li>",
        "<li>Tarasov, A., Vilella, A. J., Cuppen, E., Nijman, I. J., & Prins, P. (2015). Sambamba: Fast processing of NGS alignment formats. Bioinformatics, 31(12), 2032–2034. https://doi.org/10.1093/bioinformatics/btv098</li>",
        "<li>Eisfeldt, J., Vezzi, F., Olason, P., Nilsson, D., & Lindstrand, A. (2017). TIDDIT, an efficient and comprehensive structural variant caller for massive parallel sequencing data. F1000Research, 6, 664. https://doi.org/10.12688/f1000research.11168.2</li>",
        "<li>Kent, W. J., Zweig, A. S., Barber, G., Hinrichs, A. S., & Karolchik, D. (2010). BigWig and BigBed: Enabling browsing of large distributed datasets. Bioinformatics, 26(17), 2204–2207. https://doi.org/10.1093/bioinformatics/btq351</li>",
        (params.verifybamid_svd_bed && params.verifybamid_svd_mu && params.verifybamid_svd_ud) ? "<li>Zhang, F., Flickinger, M., Taliun, S. A. G., Consortium, I. P. G., Abecasis, G. R., Scott, L. J., McCaroll, S. A., Pato, C. N., Boehnke, M., & Kang, H. M. (2020). Ancestry-agnostic estimation of DNA sample contamination from sequence reads. Genome Research, 30(2), 185–194. https://doi.org/10.1101/gr.246934.118</li>" : "",
        "<li>Pedersen, B. S., & Quinlan, A. R. (2018). Mosdepth: Quick coverage calculation for genomes and exomes. Bioinformatics, 34(5), 867–868. https://doi.org/10.1093/bioinformatics/btx699</li>"
    ]
    preprocessing_text = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</li>",
        (params.skip_tools && params.skip_tools.split(',').contains('fastp')) ? "" : "<li>Chen, S. (2023). Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta, 2(2), e107. https://doi.org/10.1002/imt2.107</li>",
        hasSpringInput() ? "<li>Chandak, S., Tatwawadi, K., Ochoa, I., Hernaez, M., & Weissman, T. (2019). SPRING: A next-generation compressor for FASTQ data. Bioinformatics, 35(15), 2674–2676. https://doi.org/10.1093/bioinformatics/bty1015</li>" : ""
    ]

    other_citation_text = [
        "<li>Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008</li>",
        "<li>McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., & DePristo, M. A. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297–1303. https://doi.org/10.1101/gr.107524.110</li>",
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354</li>",
        (params.skip_tools && params.skip_tools.split(',').contains('peddy')) ? "" : "<li>Pedersen, B. S., & Quinlan, A. R. (2017). Who’s Who? Detecting and Resolving Sample Anomalies in Human DNA Sequencing Studies with Peddy. The American Journal of Human Genetics, 100(3), 406–413. https://doi.org/10.1016/j.ajhg.2017.01.017</li>",
        params.run_rtgvcfeval ? "<li>Cleary, J. G., Braithwaite, R., Gaastra, K., Hilbush, B. S., Inglis, S., Irvine, S. A., Jackson, A., Littin, R., Rathod, M., Ware, D., Zook, J. M., Trigg, L., & Vega, F. M. D. L. (2015). Comparing Variant Call Files for Performance Benchmarking of Next-Generation Sequencing Variant Calling Pipelines (p. 023754). bioRxiv. https://doi.org/10.1101/023754</li>" : "",
        "<li>Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352</li>",
        (!(params.skip_tools && params.skip_tools.split(',').contains('smncopynumbercaller')) && params.analysis_type.equals("wgs")) ? "<li>Chen, X., Sanchis-Juan, A., French, C. E., Connell, A. J., Delon, I., Kingsbury, Z., Chawla, A., Halpern, A. L., Taft, R. J., Bentley, D. R., Butchbach, M. E. R., Raymond, F. L., & Eberle, M. A. (2020). Spinal muscular atrophy diagnosis and carrier screening from genome sequencing data. Genetics in Medicine, 22(5), 945–953. https://doi.org/10.1038/s41436-020-0754-0</li>" : "",
        "<li>Li, H. (2011). Tabix: Fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics, 27(5), 718–719. https://doi.org/10.1093/bioinformatics/btq671</li>",
        "<li>Quinlan, AR., Hall IM. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinfomatics, 26(6), 841-842. https://doi.org/10.1093/bioinformatics/btq033</li>"
    ]

    def concat_text = align_text +
                        repeats_text         +
                        snv_annotation_text  +
                        snv_calls_text       +
                        sv_annotation_text   +
                        sv_calls_text        +
                        mt_annotation_text   +
                        qc_bam_text          +
                        me_calls_text        +
                        me_annotation_text   +
                        preprocessing_text   +
                        other_citation_text

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
