/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowRaredisease.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CHECK MANDATORY PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def mandatoryParams = [
    "aligner",
    "analysis_type",
    "fasta",
    "input",
    "intervals_wgs",
    "intervals_y",
    "platform",
    "variant_catalog",
    "variant_caller"
]

if (!params.skip_snv_annotation) {
    mandatoryParams += ["genome", "vcfanno_resources", "vcfanno_toml", "vep_cache", "vep_cache_version",
    "gnomad_af", "score_config_snv"]
}

if (!params.skip_sv_annotation) {
    mandatoryParams += ["genome", "svdb_query_dbs", "vep_cache", "vep_cache_version", "score_config_sv"]
}

if (!params.skip_mt_annotation) {
    mandatoryParams += ["genome", "mito_name", "vcfanno_resources", "vcfanno_toml", "vep_cache_version", "vep_cache"]
}

if (params.analysis_type.equals("wes")) {
    mandatoryParams += ["target_bed"]
}

if (params.variant_caller.equals("sentieon")) {
    mandatoryParams += ["ml_model"]
}

if (!params.skip_germlinecnvcaller) {
    mandatoryParams += ["ploidy_model", "gcnvcaller_model"]
}

if (!params.skip_vep_filter) {
    mandatoryParams += ["vep_filters"]
}

def missingParamsCount = 0
for (param in mandatoryParams.unique()) {
    if (params[param] == null) {
        println("params." + param + " not set.")
        missingParamsCount += 1
    }
}

if (missingParamsCount>0) {
    error("\nSet missing parameters and restart the run. For more information please check usage documentation on github.")
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config              = params.multiqc_config              ? Channel.fromPath( params.multiqc_config )  : Channel.empty()
ch_multiqc_logo                       = params.multiqc_logo                ? Channel.fromPath( params.multiqc_logo )    : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description )  : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//
// MODULE: Installed directly from nf-core/modules
//

include { CUSTOM_DUMPSOFTWAREVERSIONS           } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTQC                                } from '../modules/nf-core/fastqc/main'
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'
include { SMNCOPYNUMBERCALLER                   } from '../modules/nf-core/smncopynumbercaller/main'
include { ENSEMBLVEP_FILTERVEP as FILTERVEP_MT } from '../modules/nf-core/ensemblvep/filtervep'
include { ENSEMBLVEP_FILTERVEP as FILTERVEP_SNV } from '../modules/nf-core/ensemblvep/filtervep'
include { ENSEMBLVEP_FILTERVEP as FILTERVEP_SV  } from '../modules/nf-core/ensemblvep/filtervep'
include { TABIX_BGZIPTABIX as BGZIPTABIX_MT    } from '../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as BGZIPTABIX_SNV    } from '../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as BGZIPTABIX_SV     } from '../modules/nf-core/tabix/bgziptabix'

//
// SUBWORKFLOWS
//

include { ALIGN                                 } from '../subworkflows/local/align'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_MT    } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SNV   } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SV    } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_GENOME_SNVS                  } from '../subworkflows/local/annotate_genome_snvs'
include { ANNOTATE_MT_SNVS                      } from '../subworkflows/local/annotate_mt_snvs'
include { ANNOTATE_STRUCTURAL_VARIANTS          } from '../subworkflows/local/annotate_structural_variants'
include { CALL_REPEAT_EXPANSIONS                } from '../subworkflows/local/call_repeat_expansions'
include { CALL_SNV                              } from '../subworkflows/local/call_snv'
include { CALL_STRUCTURAL_VARIANTS              } from '../subworkflows/local/call_structural_variants'
include { GENS                                  } from '../subworkflows/local/gens'
include { PREPARE_REFERENCES                    } from '../subworkflows/local/prepare_references'
include { QC_BAM                                } from '../subworkflows/local/qc_bam'
include { RANK_VARIANTS as RANK_VARIANTS_MT     } from '../subworkflows/local/rank_variants'
include { RANK_VARIANTS as RANK_VARIANTS_SNV    } from '../subworkflows/local/rank_variants'
include { RANK_VARIANTS as RANK_VARIANTS_SV     } from '../subworkflows/local/rank_variants'
include { SCATTER_GENOME                        } from '../subworkflows/local/scatter_genome'
include { PEDDY_CHECK                           } from '../subworkflows/local/peddy_check'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow RAREDISEASE {

    ch_versions = Channel.empty()

    // Initialize read, sample, and case_info channels
    ch_input = Channel.fromPath(params.input)
    Channel.fromSamplesheet("input")
        .tap { ch_original_input }
        .map { meta, fastq1, fastq2 -> meta.id }
        .reduce([:]) { counts, sample -> //get counts of each sample in the samplesheet - for groupTuple
            counts[sample] = (counts[sample] ?: 0) + 1
            counts
        }
        .combine( ch_original_input )
        .map { counts, meta, fastq1, fastq2 ->
            new_meta = meta + [num_lanes:counts[meta.id],
                        read_group:"\'@RG\\tID:"+ fastq1.toString().split('/')[-1] + "\\tPL:ILLUMINA\\tSM:"+meta.id+"\'"]
            if (!fastq2) {
                return [ new_meta + [ single_end:true ], [ fastq1 ] ]
            } else {
                return [ new_meta + [ single_end:false ], [ fastq1, fastq2 ] ]
            }
        }
        .tap{ ch_input_counts }
        .map { meta, fastqs -> fastqs }
        .reduce([:]) { counts, fastqs -> //get line number for each row to construct unique sample ids
            counts[fastqs] = counts.size() + 1
            return counts
        }
        .combine( ch_input_counts )
        .map { lineno, meta, fastqs -> //append line number to sampleid
            new_meta = meta + [id:meta.id+"_T"+lineno[fastqs]]
            return [ new_meta, fastqs ]
        }
        .set { ch_reads }

    ch_samples   = ch_reads.map { meta, fastqs -> meta}
    ch_pedfile   = ch_samples.toList().map { makePed(it) }
    ch_case_info = ch_samples.toList().map { create_case_channel(it) }

    // Initialize file channels for PREPARE_REFERENCES subworkflow
    ch_genome_fasta             = Channel.fromPath(params.fasta).map { it -> [[id:it[0].simpleName], it] }.collect()
    ch_genome_fai               = params.fai            ? Channel.fromPath(params.fai).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                        : Channel.empty()
    ch_gnomad_af_tab            = params.gnomad_af      ? Channel.fromPath(params.gnomad_af).map{ it -> [[id:it[0].simpleName], it] }.collect()
                                                        : Channel.value([[],[]])
    ch_dbsnp                    = params.known_dbsnp    ? Channel.fromPath(params.known_dbsnp).map{ it -> [[id:it[0].simpleName], it] }.collect()
                                                        : Channel.value([[],[]])
    ch_mt_fasta                 = params.mt_fasta       ? Channel.fromPath(params.mt_fasta).map { it -> [[id:it[0].simpleName], it] }.collect()
                                                        : Channel.empty()
    ch_target_bed_unprocessed   = params.target_bed     ? Channel.fromPath(params.target_bed).map{ it -> [[id:it[0].simpleName], it] }.collect()
                                                        : Channel.value([[],[]])
    ch_vep_cache_unprocessed    = params.vep_cache      ? Channel.fromPath(params.vep_cache).map { it -> [[id:'vep_cache'], it] }.collect()
                                                        : Channel.value([[],[]])

    // Prepare references and indices.
    PREPARE_REFERENCES (
        ch_genome_fasta,
        ch_genome_fai,
        ch_mt_fasta,
        ch_gnomad_af_tab,
        ch_dbsnp,
        ch_target_bed_unprocessed,
        ch_vep_cache_unprocessed
    )
    .set { ch_references }

    // Gather built indices or get them from the params
    ch_bait_intervals           = ch_references.bait_intervals
    ch_cadd_header              = Channel.fromPath("$projectDir/assets/cadd_to_vcf_header_-1.0-.txt", checkIfExists: true).collect()
    ch_cadd_resources           = params.cadd_resources                    ? Channel.fromPath(params.cadd_resources).collect()
                                                                           : Channel.value([])
    ch_call_interval            = params.call_interval                     ? Channel.fromPath(params.call_interval).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                           : Channel.value([[:],[]])
    ch_dbsnp_tbi                = params.known_dbsnp_tbi                   ? Channel.fromPath(params.known_dbsnp_tbi).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                           : ch_references.known_dbsnp_tbi.ifEmpty([[],[]])
    ch_gcnvcaller_model         = params.gcnvcaller_model                  ? Channel.fromPath(params.gcnvcaller_model).splitCsv ( header:true )
                                                                            .map { row ->
                                                                                return [[id:file(row.models).simpleName], row.models]
                                                                            }
                                                                           : Channel.empty()
    ch_genome_bwaindex          = params.bwa                               ? Channel.fromPath(params.bwa).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                           : ch_references.genome_bwa_index
    ch_genome_bwamem2index      = params.bwamem2                           ? Channel.fromPath(params.bwamem2).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                           : ch_references.genome_bwamem2_index
    ch_genome_chrsizes          = ch_references.genome_chrom_sizes
    ch_genome_fai               = ch_references.genome_fai
    ch_genome_dictionary        = params.sequence_dictionary               ? Channel.fromPath(params.sequence_dictionary).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                           : ch_references.genome_dict
    ch_gnomad_afidx             = params.gnomad_af_idx                     ? Channel.fromPath(params.gnomad_af_idx).collect()
                                                                           : ch_references.gnomad_af_idx
    ch_gnomad_af                = params.gnomad_af                         ? ch_gnomad_af_tab.join(ch_gnomad_afidx).map {meta, tab, idx -> [tab,idx]}.collect()
                                                                           : Channel.empty()
    ch_intervals_wgs            = params.intervals_wgs                     ? Channel.fromPath(params.intervals_wgs).collect()
                                                                           : Channel.empty()
    ch_intervals_y              = params.intervals_y                       ? Channel.fromPath(params.intervals_y).collect()
                                                                           : Channel.empty()
    ch_ml_model                 = params.variant_caller.equals("sentieon") ? Channel.fromPath(params.ml_model).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                           : Channel.value([[:],[]])
    ch_mt_intervals             = ch_references.mt_intervals
    ch_mtshift_backchain        = ch_references.mtshift_backchain
    ch_mtshift_bwaindex         = ch_references.mtshift_bwa_index
    ch_mtshift_bwamem2index     = ch_references.mtshift_bwamem2_index
    ch_mtshift_dictionary       = ch_references.mtshift_dict
    ch_mtshift_fai              = ch_references.mtshift_fai
    ch_mtshift_fasta            = ch_references.mtshift_fasta
    ch_mtshift_intervals        = ch_references.mtshift_intervals
    ch_ploidy_model             = params.ploidy_model                      ? Channel.fromPath(params.ploidy_model).map{ it -> [[id:it[0].simpleName], it] }.collect()
                                                                           : Channel.empty()
    ch_readcount_intervals      = params.readcount_intervals               ? Channel.fromPath(params.readcount_intervals).collect()
                                                                           : ( ch_references.readcount_intervals      ?: Channel.empty() )
    ch_reduced_penetrance       = params.reduced_penetrance                ? Channel.fromPath(params.reduced_penetrance).collect()
                                                                           : Channel.value([])
    ch_score_config_mt          = params.score_config_mt                   ? Channel.fromPath(params.score_config_mt).collect()
                                                                           : Channel.value([])
    ch_score_config_snv         = params.score_config_snv                  ? Channel.fromPath(params.score_config_snv).collect()
                                                                           : Channel.value([])
    ch_score_config_sv          = params.score_config_sv                   ? Channel.fromPath(params.score_config_sv).collect()
                                                                           : Channel.value([])
    ch_target_bed               = ch_references.target_bed
    ch_target_intervals         = ch_references.target_intervals
    ch_variant_catalog          = params.variant_catalog                   ? Channel.fromPath(params.variant_catalog).map { it -> [[id:it[0].simpleName],it]}.collect()
                                                                           : Channel.value([[],[]])
    ch_variant_consequences     = Channel.fromPath("$projectDir/assets/variant_consequences_v1.txt", checkIfExists: true).collect()
    ch_vcfanno_resources        = params.vcfanno_resources                 ? Channel.fromPath(params.vcfanno_resources).splitText().map{it -> it.trim()}.collect()
                                                                           : Channel.value([])
    ch_vcfanno_lua              = params.vcfanno_lua                       ? Channel.fromPath(params.vcfanno_lua).collect()
                                                                           : Channel.value([])
    ch_vcfanno_toml             = params.vcfanno_toml                      ? Channel.fromPath(params.vcfanno_toml).collect()
                                                                           : Channel.value([])
    ch_vep_cache                = ( params.vep_cache && params.vep_cache.endsWith("tar.gz") )  ? ch_references.vep_resources
                                                                           : ( params.vep_cache    ? Channel.fromPath(params.vep_cache).collect() : Channel.value([]) )
    ch_vep_filters              = params.vep_filters                       ? Channel.fromPath(params.vep_filters).collect()
                                                                           : Channel.value([])
    ch_versions                 = ch_versions.mix(ch_references.versions)


    // SV caller priority
    if (params.skip_germlinecnvcaller) {
        ch_svcaller_priority = Channel.value(["tiddit", "manta"])
    } else {
        ch_svcaller_priority = Channel.value(["tiddit", "manta", "gcnvcaller"])
    }

    // Input QC
    if (!params.skip_fastqc) {
        FASTQC (ch_reads)
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }
    // CREATE CHROMOSOME BED AND INTERVALS
    SCATTER_GENOME (
        ch_genome_dictionary,
        ch_genome_fai,
        ch_genome_fasta
    )
    .set { ch_scatter }

    ch_scatter_split_intervals  = ch_scatter.split_intervals  ?: Channel.empty()

    //
    // ALIGNING READS, FETCH STATS, AND MERGE.
    //
    ALIGN (
        ch_reads,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_bwaindex,
        ch_genome_bwamem2index,
        ch_genome_dictionary,
        ch_mtshift_bwaindex,
        ch_mtshift_bwamem2index,
        ch_mtshift_fasta,
        ch_mtshift_dictionary,
        ch_mtshift_fai,
        params.platform
    )
    .set { ch_mapped }
    ch_versions   = ch_versions.mix(ALIGN.out.versions)

    //
    // BAM QUALITY CHECK
    //
    QC_BAM (
        ch_mapped.genome_marked_bam,
        ch_mapped.genome_marked_bai,
        ch_mapped.genome_bam_bai,
        ch_genome_fasta,
        ch_genome_fai,
        ch_bait_intervals,
        ch_target_intervals,
        ch_genome_chrsizes,
        ch_intervals_wgs,
        ch_intervals_y
    )
    ch_versions = ch_versions.mix(QC_BAM.out.versions)

    //
    // EXPANSIONHUNTER AND STRANGER
    //
    CALL_REPEAT_EXPANSIONS (
        ch_mapped.genome_bam_bai,
        ch_variant_catalog,
        ch_case_info,
        ch_genome_fasta,
        ch_genome_fai
    )
    ch_versions = ch_versions.mix(CALL_REPEAT_EXPANSIONS.out.versions)

    //
    // SNV CALLING
    //
    CALL_SNV (
        ch_mapped.genome_bam_bai,
        ch_mapped.mt_bam_bai,
        ch_mapped.mtshift_bam_bai,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_dictionary,
        ch_mt_intervals,
        ch_mtshift_fasta,
        ch_mtshift_fai,
        ch_mtshift_dictionary,
        ch_mtshift_intervals,
        ch_mtshift_backchain,
        ch_dbsnp,
        ch_dbsnp_tbi,
        ch_call_interval,
        ch_ml_model,
        ch_case_info,
        Channel.value(params.sentieon_dnascope_pcr_indel_model)
    )
    ch_versions = ch_versions.mix(CALL_SNV.out.versions)

    //
    // SV CALLING
    //
    CALL_STRUCTURAL_VARIANTS (
        ch_mapped.genome_marked_bam,
        ch_mapped.genome_marked_bai,
        ch_mapped.genome_bam_bai,
        ch_mapped.mt_bam_bai,
        ch_mapped.mtshift_bam_bai,
        ch_genome_bwaindex,
        ch_genome_fasta,
        ch_genome_fai,
        ch_mtshift_fasta,
        ch_case_info,
        ch_target_bed,
        ch_genome_dictionary,
        ch_svcaller_priority,
        ch_readcount_intervals,
        ch_ploidy_model,
        ch_gcnvcaller_model
    )
    ch_versions = ch_versions.mix(CALL_STRUCTURAL_VARIANTS.out.versions)

    //
    // ANNOTATE STRUCTURAL VARIANTS
    //
    if (!params.skip_sv_annotation) {
        ANNOTATE_STRUCTURAL_VARIANTS (
            CALL_STRUCTURAL_VARIANTS.out.vcf,
            params.svdb_query_dbs,
            params.genome,
            params.vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta,
            ch_genome_dictionary
        ).set {ch_sv_annotate}
        ch_versions = ch_versions.mix(ch_sv_annotate.versions)

        ANN_CSQ_PLI_SV (
            ch_sv_annotate.vcf_ann,
            ch_variant_consequences
        )
        ch_versions = ch_versions.mix(ANN_CSQ_PLI_SV.out.versions)

        RANK_VARIANTS_SV (
            ANN_CSQ_PLI_SV.out.vcf_ann,
            ch_pedfile,
            ch_reduced_penetrance,
            ch_score_config_sv
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_SV.out.versions)

        FILTERVEP_SV(
            RANK_VARIANTS_SV.out.vcf,
            ch_vep_filters
        )
        ch_versions = ch_versions.mix(FILTERVEP_SV.out.versions)

        BGZIPTABIX_SV(FILTERVEP_SV.out.output)
        ch_versions = ch_versions.mix(BGZIPTABIX_SV.out.versions)

    }

    //
    // ANNOTATE GENOME SNVs
    //
    if (!params.skip_snv_annotation) {

        ANNOTATE_GENOME_SNVS (
            CALL_SNV.out.genome_vcf_tabix,
            params.analysis_type,
            ch_cadd_header,
            ch_cadd_resources,
            ch_vcfanno_resources,
            ch_vcfanno_lua,
            ch_vcfanno_toml,
            params.genome,
            params.vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta,
            ch_gnomad_af,
            ch_scatter_split_intervals
        ).set {ch_snv_annotate}
        ch_versions = ch_versions.mix(ch_snv_annotate.versions)

        ch_snv_annotate = ANNOTATE_GENOME_SNVS.out.vcf_ann

        ANN_CSQ_PLI_SNV (
            ch_snv_annotate,
            ch_variant_consequences
        )
        ch_versions = ch_versions.mix(ANN_CSQ_PLI_SNV.out.versions)

        RANK_VARIANTS_SNV (
            ANN_CSQ_PLI_SNV.out.vcf_ann,
            ch_pedfile,
            ch_reduced_penetrance,
            ch_score_config_snv
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_SNV.out.versions)

        FILTERVEP_SNV(
            RANK_VARIANTS_SNV.out.vcf,
            ch_vep_filters
        )
        ch_versions = ch_versions.mix(FILTERVEP_SNV.out.versions)

        BGZIPTABIX_SNV(FILTERVEP_SNV.out.output)
        ch_versions = ch_versions.mix(BGZIPTABIX_SNV.out.versions)

    }

    //
    // ANNOTATE MT SNVs
    //
    if (!params.skip_mt_annotation) {

        ANNOTATE_MT_SNVS (
            CALL_SNV.out.mt_vcf,
            CALL_SNV.out.mt_tabix,
            ch_cadd_header,
            ch_cadd_resources,
            ch_genome_fasta,
            ch_vcfanno_resources,
            ch_vcfanno_toml,
            params.genome,
            params.vep_cache_version,
            ch_vep_cache,
        ).set {ch_mt_annotate}
        ch_versions = ch_versions.mix(ch_mt_annotate.versions)

        ANN_CSQ_PLI_MT (
            ch_mt_annotate.vcf_ann,
            ch_variant_consequences
        )
        ch_versions = ch_versions.mix(ANN_CSQ_PLI_MT.out.versions)

        RANK_VARIANTS_MT (
            ANN_CSQ_PLI_MT.out.vcf_ann,
            ch_pedfile,
            ch_reduced_penetrance,
            ch_score_config_mt
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_MT.out.versions)

        FILTERVEP_MT(
            RANK_VARIANTS_MT.out.vcf,
            ch_vep_filters
        )
        ch_versions = ch_versions.mix(FILTERVEP_MT.out.versions)

        BGZIPTABIX_MT(FILTERVEP_MT.out.output)
        ch_versions = ch_versions.mix(BGZIPTABIX_MT.out.versions)

    }

    // STEP 1.7: SMNCOPYNUMBERCALLER
    ch_mapped.genome_bam_bai
        .collect{it[1]}
        .toList()
        .set { ch_bam_list }

    ch_mapped.genome_bam_bai
        .collect{it[2]}
        .toList()
        .set { ch_bai_list }

    ch_case_info
        .combine(ch_bam_list)
        .combine(ch_bai_list)
        .set { ch_bams_bais }

    SMNCOPYNUMBERCALLER (
        ch_bams_bais
    )
    ch_versions = ch_versions.mix(SMNCOPYNUMBERCALLER.out.versions)

    // ped correspondence, sex check, ancestry check
    PEDDY_CHECK (
        CALL_SNV.out.genome_vcf.join(CALL_SNV.out.genome_tabix, failOnMismatch:true, failOnDuplicate:true),
        ch_pedfile
    )
    ch_versions = ch_versions.mix(PEDDY_CHECK.out.versions)

    // GENS
    if (params.gens_switch) {
        GENS (
            ch_mapped.genome_bam_bai,
            CALL_SNV.out.vcf,
            ch_genome_fasta,
            ch_genome_fai,
            file(params.gens_interval_list),
            file(params.gens_pon),
            file(params.gens_gnomad_pos),
            ch_case_info,
            ch_genome_dictionary
        )
        ch_versions = ch_versions.mix(GENS.out.versions)
    }

    //
    // MODULE: Pipeline reporting
    //

    // The template v2.7.1 template update introduced: ch_versions.unique{ it.text }.collectFile(name: 'collated_versions.yml')
    // This caused the pipeline to stall
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowRaredisease.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowRaredisease.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    if (!params.skip_fastqc) {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    }
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.multiple_metrics.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.hs_metrics.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.qualimap_results.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.global_dist.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.cov.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PEDDY_CHECK.out.ped.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PEDDY_CHECK.out.csv.map{it[1]}.collect().ifEmpty([]))


    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
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
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def makePed(samples) {

    def case_name  = samples[0].case_id
    def outfile  = file("${params.outdir}/pipeline_info/${case_name}" + '.ped')
    outfile.text = ['#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype'].join('\t')
    def samples_list = []
    for(int i = 0; i<samples.size(); i++) {
        sample_name        =  samples[i].sample
        if (!samples_list.contains(sample_name)) {
            outfile.append('\n' + [samples[i].case_id, sample_name, samples[i].paternal, samples[i].maternal, samples[i].sex, samples[i].phenotype].join('\t'));
            samples_list.add(sample_name)
        }
    }
    return outfile
}

// Function to get a list of metadata (e.g. case id) for the case [ meta ]
def create_case_channel(List rows) {
    def case_info    = [:]
    def probands     = []
    def upd_children = []
    def father       = ""
    def mother       = ""

    for (item in rows) {
        if (item.phenotype == "2") {
            probands.add(item.sample)
        }
        if ( (item.paternal!="0") && (item.paternal!="") && (item.maternal!="0") && (item.maternal!="") ) {
            upd_children.add(item.sample)
        }
        if ( (item.paternal!="0") && (item.paternal!="") ) {
            father = item.paternal
        }
        if ( (item.maternal!="0") && (item.maternal!="") ) {
            mother = item.maternal
        }
    }

    case_info.father       = father
    case_info.mother       = mother
    case_info.probands     = probands.unique()
    case_info.upd_children = upd_children.unique()
    case_info.id           = rows[0].case_id

    return case_info
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
