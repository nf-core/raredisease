/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRaredisease.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.bwamem2_index,
    params.fasta,
    params.fasta_fai,
    params.gnomad,
    params.input,
    params.multiqc_config,
    params.vcfanno_resources
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK    } from '../subworkflows/local/input_check'
include { CHECK_VCF      } from '../subworkflows/local/prepare_vcf'
include { CHECK_BED      } from '../subworkflows/local/prepare_bed'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

//
// SUBWORKFLOW: Consists entirely of nf-core/modules
//

include { ALIGN                       } from  '../subworkflows/nf-core/align'
include { CALL_REPEAT_EXPANSIONS      } from '../subworkflows/nf-core/call_repeat_expansions'
include { CALL_SNV_DEEPVARIANT        } from '../subworkflows/nf-core/call_snv_deepvariant'
include { QC_BAM                      } from '../subworkflows/nf-core/qc_bam'
include { ANNOTATE_VCFANNO            } from '../subworkflows/nf-core/annotate_vcfanno'
include { CALL_STRUCTURAL_VARIANTS    } from '../subworkflows/nf-core/call_structural_variants'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RAREDISEASE {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // STEP 0: QUALITY CHECK.
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // STEP 0: PREPARE GENOME REFERENCES AND INDICES.
    PREPARE_GENOME (
        params.fasta,
        params.variant_catalog
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    ch_gnomad = Channel.empty()
    if (params.gnomad) {
        CHECK_VCF(
            params.gnomad,
            PREPARE_GENOME.out.fasta
        )
        .set { ch_gnomad }
    }

    ch_target_bed = Channel.empty()
    if (params.target_bed) {
        CHECK_BED(
            params.target_bed,
            PREPARE_GENOME.out.sequence_dict
        )
        .set { ch_target_bed }
    }

    // STEP 1: ALIGNING READS, FETCH STATS, AND MERGE.
    ALIGN (
        params.aligner,
        INPUT_CHECK.out.reads,
        PREPARE_GENOME.out.bwamem2_index
    )
    ch_marked_bam = ALIGN.out.marked_bam
    ch_marked_bai = ALIGN.out.marked_bai
    ch_bam_bai    = ch_marked_bam.join(ch_marked_bai, by: [0])

    ch_versions   = ch_versions.mix(ALIGN.out.versions)

    // STEP 1.5: BAM QUALITY CHECK
    QC_BAM (
        ch_marked_bam,
        ch_marked_bai,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        CHECK_BED.out.bait_intervals,
        CHECK_BED.out.target_intervals,
        PREPARE_GENOME.out.chrom_sizes
    )
    ch_versions = ch_versions.mix(QC_BAM.out.versions.ifEmpty(null))

    // STEP 1.6: EXPANSIONHUNTER
    CALL_REPEAT_EXPANSIONS (
        ch_bam_bai,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.variant_catalog
    )
    ch_versions = ch_versions.mix(CALL_REPEAT_EXPANSIONS.out.versions.ifEmpty(null))

    // STEP 2: VARIANT CALLING
    // TODO: There should be a conditional to execute certain variant callers (e.g. sentieon, gatk, deepvariant) defined by the user and we need to think of a default caller.
    CALL_SNV_DEEPVARIANT (
        ch_bam_bai,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        INPUT_CHECK.out.case_info
    )
    ch_versions = ch_versions.mix(CALL_SNV_DEEPVARIANT.out.versions)

    CALL_STRUCTURAL_VARIANTS (
        ch_marked_bam,
        ch_marked_bai,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        INPUT_CHECK.out.case_info,
        ch_target_bed.bed
    )
    ch_versions = ch_versions.mix(CALL_STRUCTURAL_VARIANTS.out.versions)

    // STEP 3: VARIANT ANNOTATION
    ch_dv_vcf = CALL_SNV_DEEPVARIANT.out.vcf.join(CALL_SNV_DEEPVARIANT.out.tabix, by: [0])

    ANNOTATE_VCFANNO (
        params.vcfanno_toml,
        ch_dv_vcf,
        PREPARE_GENOME.out.vcfanno_resources
    )
    ch_versions = ch_versions.mix(ANNOTATE_VCFANNO.out.versions)

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowRaredisease.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
