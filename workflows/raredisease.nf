/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRaredisease.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.multiqc_config, params.fasta,
    params.bwamem2_index, params.fasta_fai
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

// Switch for saving references
if (!params.save_reference) {
    modules['bwa_mem2_index'].publish_files = false
    modules['samtools_faidx'].publish_files = false
}

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )

//
// SUBWORKFLOW: Consists entirely of nf-core/modules
//

include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome' addParams(
    bwamem2_idx_options: modules['bwa_mem2_index'],
    samtools_faidx_options: modules['samtools_faidx']
)

include { ALIGN_BWAMEM2 } from  '../subworkflows/nf-core/align_bwamem2' addParams(
    bwamem2_idx_options: modules['bwa_mem2_index'],
    bwamem2_mem_options: modules['bwa_mem2_mem'],
    samtools_idx_options: modules['samtools_index'],
    samtools_sort_options: modules['samtools_sort'],
    samtools_stats_options: modules['samtools_stats'],
    samtools_merge_options: modules['samtools_merge'],
    markduplicates_options: modules['picard_markduplicates']
)

include { QC_BAM } from '../subworkflows/nf-core/qc_bam' addParams (
    picard_collectmultiplemetrics_options: modules['picard_collectmultiplemetrics'],
    qualimap_bamqc_options: modules['qualimap_bamqc']
)


//
// SUBWORKFLOW: Consists of mix/local modules
//

include { DEEPVARIANT_CALLER } from '../subworkflows/local/deepvariant_caller' addParams( deepvariant_options: modules['deepvariant'],
                                                                                        glnexus_options: modules['glnexus'],
                                                                                        rm_duplicates_options: modules['bcftools_norm_rm_duplicates'],
                                                                                        split_multiallelics_options: modules['bcftools_norm_split_multiallelics'],
                                                                                        tabix_options: modules['tabix'] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RAREDISEASE {

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )

    // STEP 0: QUALITY CHECK.
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.ifEmpty(null))

    // STEP 0: PREPARE GENOME REFERENCES AND INDICES.
    PREPARE_GENOME ( params.fasta )

    // STEP 1: ALIGNING READS, FETCH STATS, AND MERGE.
    if (params.aligner == 'bwamem2') {
        ALIGN_BWAMEM2 (
            INPUT_CHECK.out.reads,
            PREPARE_GENOME.out.bwamem2_index
        )

        ch_marked_bam = ALIGN_BWAMEM2.out.marked_bam
        ch_marked_bai = ALIGN_BWAMEM2.out.marked_bai

        ch_software_versions = ch_software_versions.mix(ALIGN_BWAMEM2.out.bwamem2_version.ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ALIGN_BWAMEM2.out.picard_version.ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(ALIGN_BWAMEM2.out.samtools_version.ifEmpty(null))
    }

    // STEP 1.5: BAM QUALITY CHECK
    gff = []
    use_gff = false
    QC_BAM (
        ch_marked_bam,
        PREPARE_GENOME.out.fasta,
        gff,
        use_gff
    )

    // STEP 2: VARIANT CALLING
    // TODO: There should be a conditional to execute certain variant callers (e.g. sentieon, gatk, deepvariant) defined by the user and we need to think of a default caller.
    DEEPVARIANT_CALLER (
                        ch_marked_bam.join(ch_marked_bai, by: [0]),
                        PREPARE_GENOME.out.fasta,
                        PREPARE_GENOME.out.fai,
                        INPUT_CHECK.out.ch_case_info
                        )
    ch_software_versions = ch_software_versions.mix(DEEPVARIANT_CALLER.out.deepvariant_version.ifEmpty(null))

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowRaredisease.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files    = Channel.empty()
    ch_multiqc_files    = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files    = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files    = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files    = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files    = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
