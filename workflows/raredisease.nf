/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_raredisease_pipeline'

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
def missingParamsCount = 0

if (params.run_rtgvcfeval) {
    mandatoryParams += ["rtg_truthvcfs"]
}

if (!params.skip_snv_annotation) {
    mandatoryParams += ["genome", "vcfanno_resources", "vcfanno_toml", "vep_cache", "vep_cache_version",
    "gnomad_af", "score_config_snv", "variant_consequences_snv"]
}

if (!params.skip_sv_annotation) {
    mandatoryParams += ["genome", "vep_cache", "vep_cache_version", "score_config_sv", "variant_consequences_sv"]
    if (!params.svdb_query_bedpedbs && !params.svdb_query_dbs) {
        println("params.svdb_query_bedpedbs or params.svdb_query_dbs should be set.")
        missingParamsCount += 1
    }
}

if (!params.skip_mt_annotation) {
    mandatoryParams += ["genome", "mito_name", "vcfanno_resources", "vcfanno_toml", "vep_cache_version", "vep_cache", "variant_consequences_snv"]
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
    if (!params.vep_filters && !params.vep_filters_scout_fmt) {
        println("params.vep_filters or params.vep_filters_scout_fmt should be set.")
        missingParamsCount += 1
    } else if (params.vep_filters && params.vep_filters_scout_fmt) {
        println("Either params.vep_filters or params.vep_filters_scout_fmt should be set.")
        missingParamsCount += 1
    }
}

if (!params.skip_me_annotation) {
    mandatoryParams += ["mobile_element_svdb_annotations", "variant_consequences_snv"]
}

if (!params.skip_gens) {
    mandatoryParams += ["gens_gnomad_pos", "gens_interval_list", "gens_pon_female", "gens_pon_male"]
}

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
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//
// MODULE: Installed directly from nf-core/modules
//

include { FASTQC              } from '../modules/nf-core/fastqc/main'
include { MULTIQC             } from '../modules/nf-core/multiqc/main'
include { PEDDY               } from '../modules/nf-core/peddy/main'
include { SMNCOPYNUMBERCALLER } from '../modules/nf-core/smncopynumbercaller/main'

//
// MODULE: Local modules
//

include { RENAME_ALIGN_FILES as RENAME_BAM_FOR_SMNCALLER } from '../modules/local/rename_align_files'
include { RENAME_ALIGN_FILES as RENAME_BAI_FOR_SMNCALLER } from '../modules/local/rename_align_files'

//
// SUBWORKFLOWS
//

include { ALIGN                                              } from '../subworkflows/local/align'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_MT                 } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SNV                } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SV                 } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_GENOME_SNVS                               } from '../subworkflows/local/annotate_genome_snvs'
include { ANNOTATE_MOBILE_ELEMENTS                           } from '../subworkflows/local/annotate_mobile_elements'
include { ANNOTATE_MT_SNVS                                   } from '../subworkflows/local/annotate_mt_snvs'
include { ANNOTATE_STRUCTURAL_VARIANTS                       } from '../subworkflows/local/annotate_structural_variants'
include { CALL_MOBILE_ELEMENTS                               } from '../subworkflows/local/call_mobile_elements'
include { CALL_REPEAT_EXPANSIONS                             } from '../subworkflows/local/call_repeat_expansions'
include { CALL_SNV                                           } from '../subworkflows/local/call_snv'
include { CALL_STRUCTURAL_VARIANTS                           } from '../subworkflows/local/call_structural_variants'
include { GENERATE_CLINICAL_SET as GENERATE_CLINICAL_SET_MT  } from '../subworkflows/local/generate_clinical_set'
include { GENERATE_CLINICAL_SET as GENERATE_CLINICAL_SET_SNV } from '../subworkflows/local/generate_clinical_set'
include { GENERATE_CLINICAL_SET as GENERATE_CLINICAL_SET_SV  } from '../subworkflows/local/generate_clinical_set'
include { GENERATE_CYTOSURE_FILES                            } from '../subworkflows/local/generate_cytosure_files'
include { GENS                                               } from '../subworkflows/local/gens'
include { PREPARE_REFERENCES                                 } from '../subworkflows/local/prepare_references'
include { QC_BAM                                             } from '../subworkflows/local/qc_bam'
include { RANK_VARIANTS as RANK_VARIANTS_MT                  } from '../subworkflows/local/rank_variants'
include { RANK_VARIANTS as RANK_VARIANTS_SNV                 } from '../subworkflows/local/rank_variants'
include { RANK_VARIANTS as RANK_VARIANTS_SV                  } from '../subworkflows/local/rank_variants'
include { SCATTER_GENOME                                     } from '../subworkflows/local/scatter_genome'
include { SUBSAMPLE_MT                                       } from '../subworkflows/local/subsample_mt'
include { VARIANT_EVALUATION                                 } from '../subworkflows/local/variant_evaluation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RAREDISEASE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_samples   = ch_samplesheet.map { meta, fastqs -> meta}
    ch_pedfile   = ch_samples.toList().map { file(CustomFunctions.makePed(it, params.outdir)) }
    ch_case_info = ch_samples.toList().map { CustomFunctions.createCaseChannel(it) }

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
    ch_cadd_resources           = params.cadd_resources                     ? Channel.fromPath(params.cadd_resources).collect()
                                                                            : Channel.value([])
    ch_call_interval            = params.call_interval                      ? Channel.fromPath(params.call_interval).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                            : Channel.value([[:],[]])
    ch_dbsnp_tbi                = params.known_dbsnp_tbi                    ? Channel.fromPath(params.known_dbsnp_tbi).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                            : ch_references.known_dbsnp_tbi.ifEmpty([[],[]])
    ch_foundin_header           = Channel.fromPath("$projectDir/assets/foundin.hdr", checkIfExists: true).collect()
    ch_gcnvcaller_model         = params.gcnvcaller_model                   ? Channel.fromPath(params.gcnvcaller_model).splitCsv ( header:true )
                                                                            .map { row ->
                                                                                return [[id:file(row.models).simpleName], row.models]
                                                                            }
                                                                            : Channel.empty()
    ch_genome_bwaindex          = params.bwa                                ? Channel.fromPath(params.bwa).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                            : ch_references.genome_bwa_index
    ch_genome_bwamem2index      = params.bwamem2                            ? Channel.fromPath(params.bwamem2).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                            : ch_references.genome_bwamem2_index
    ch_genome_chrsizes          = ch_references.genome_chrom_sizes
    ch_genome_fai               = ch_references.genome_fai
    ch_genome_dictionary        = params.sequence_dictionary                ? Channel.fromPath(params.sequence_dictionary).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                            : ch_references.genome_dict
    ch_gens_gnomad_pos          = params.gens_gnomad_pos                    ? Channel.fromPath(params.gens_gnomad_pos).collect()
                                                                            : Channel.empty()
    ch_gens_interval_list       = params.gens_interval_list                 ? Channel.fromPath(params.gens_interval_list).collect()
                                                                            : Channel.empty()
    ch_gens_pon_female          = params.gens_pon_female                    ? Channel.fromPath(params.gens_pon_female).map { it -> [ [id:it[0].simpleName], it ] }.collect()
                                                                            : Channel.empty()
    ch_gens_pon_male            = params.gens_pon_male                      ? Channel.fromPath(params.gens_pon_male).map { it -> [ [id:it[0].simpleName], it ] }.collect()
                                                                            : Channel.empty()
    ch_gnomad_afidx             = params.gnomad_af_idx                      ? Channel.fromPath(params.gnomad_af_idx).collect()
                                                                            : ch_references.gnomad_af_idx
    ch_gnomad_af                = params.gnomad_af                          ? ch_gnomad_af_tab.join(ch_gnomad_afidx).map {meta, tab, idx -> [tab,idx]}.collect()
                                                                            : Channel.empty()
    ch_intervals_wgs            = params.intervals_wgs                      ? Channel.fromPath(params.intervals_wgs).collect()
                                                                            : Channel.empty()
    ch_intervals_y              = params.intervals_y                        ? Channel.fromPath(params.intervals_y).collect()
                                                                            : Channel.empty()
    ch_me_references            = params.mobile_element_references          ? Channel.fromSamplesheet("mobile_element_references")
                                                                            : Channel.empty()
    ch_me_svdb_resources        = params.mobile_element_svdb_annotations    ? Channel.fromPath(params.mobile_element_svdb_annotations)
                                                                            : Channel.empty()
    ch_ml_model                 = params.variant_caller.equals("sentieon")  ? Channel.fromPath(params.ml_model).map {it -> [[id:it[0].simpleName], it]}.collect()
                                                                            : Channel.value([[:],[]])
    ch_mt_intervals             = ch_references.mt_intervals
    ch_mtshift_backchain        = ch_references.mtshift_backchain
    ch_mtshift_bwaindex         = ch_references.mtshift_bwa_index
    ch_mtshift_bwamem2index     = ch_references.mtshift_bwamem2_index
    ch_mtshift_dictionary       = ch_references.mtshift_dict
    ch_mtshift_fai              = ch_references.mtshift_fai
    ch_mtshift_fasta            = ch_references.mtshift_fasta
    ch_mtshift_intervals        = ch_references.mtshift_intervals
    ch_ploidy_model             = params.ploidy_model                       ? Channel.fromPath(params.ploidy_model).map{ it -> [[id:it[0].simpleName], it] }.collect()
                                                                            : Channel.empty()
    ch_readcount_intervals      = params.readcount_intervals                ? Channel.fromPath(params.readcount_intervals).collect()
                                                                            : ( ch_references.readcount_intervals      ?: Channel.empty() )
    ch_reduced_penetrance       = params.reduced_penetrance                 ? Channel.fromPath(params.reduced_penetrance).collect()
                                                                            : Channel.value([])
    ch_rtg_truthvcfs            = params.rtg_truthvcfs                      ? Channel.fromPath(params.rtg_truthvcfs).collect()
                                                                            : Channel.value([])
    ch_sample_id_map            = params.sample_id_map                      ? Channel.fromSamplesheet("sample_id_map")
                                                                            : Channel.empty()
    ch_score_config_mt          = params.score_config_mt                    ? Channel.fromPath(params.score_config_mt).collect()
                                                                            : Channel.value([])
    ch_score_config_snv         = params.score_config_snv                   ? Channel.fromPath(params.score_config_snv).collect()
                                                                            : Channel.value([])
    ch_score_config_sv          = params.score_config_sv                    ? Channel.fromPath(params.score_config_sv).collect()
                                                                            : Channel.value([])
    ch_sdf                      = params.sdf                                ? Channel.fromPath(params.sdf).map{it -> [[id:it[0].simpleName],it]}.collect()
                                                                            : ch_references.sdf
    ch_sv_dbs                   = params.svdb_query_dbs                     ? Channel.fromPath(params.svdb_query_dbs)
                                                                            : Channel.empty()
    ch_sv_bedpedbs              = params.svdb_query_bedpedbs                ? Channel.fromPath(params.svdb_query_bedpedbs)
                                                                            : Channel.empty()
    ch_target_bed               = ch_references.target_bed
    ch_target_intervals         = ch_references.target_intervals
    ch_variant_catalog          = params.variant_catalog                    ? Channel.fromPath(params.variant_catalog).map { it -> [[id:it[0].simpleName],it]}.collect()
                                                                            : Channel.value([[],[]])
    ch_variant_consequences_snv = params.variant_consequences_snv           ? Channel.fromPath(params.variant_consequences_snv).collect()
                                                                            : Channel.value([])
    ch_variant_consequences_sv  = params.variant_consequences_sv            ? Channel.fromPath(params.variant_consequences_sv).collect()
                                                                            : Channel.value([])
    ch_vcfanno_resources        = params.vcfanno_resources                  ? Channel.fromPath(params.vcfanno_resources).splitText().map{it -> it.trim()}.collect()
                                                                            : Channel.value([])
    ch_vcf2cytosure_blacklist   = params.vcf2cytosure_blacklist             ? Channel.fromPath(params.vcf2cytosure_blacklist).collect()
                                                                            : Channel.value([])
    ch_vcfanno_lua              = params.vcfanno_lua                        ? Channel.fromPath(params.vcfanno_lua).collect()
                                                                            : Channel.value([])
    ch_vcfanno_toml             = params.vcfanno_toml                       ? Channel.fromPath(params.vcfanno_toml).collect()
                                                                            : Channel.value([])
    ch_vep_cache                = ( params.vep_cache && params.vep_cache.endsWith("tar.gz") )  ? ch_references.vep_resources
                                                                            : ( params.vep_cache    ? Channel.fromPath(params.vep_cache).collect() : Channel.value([]) )
    ch_vep_extra_files_unsplit  = params.vep_plugin_files                   ? Channel.fromPath(params.vep_plugin_files).collect()
                                                                            : Channel.value([])
    ch_vep_filters_std_fmt      = params.vep_filters                        ? Channel.fromPath(params.vep_filters).splitCsv().collect()
                                                                            : Channel.empty()
    ch_vep_filters_scout_fmt    = params.vep_filters_scout_fmt              ? Channel.fromPath(params.vep_filters_scout_fmt).collect()
                                                                            : Channel.empty()
    ch_versions                 = ch_versions.mix(ch_references.versions)

    // SV caller priority
    if (params.skip_germlinecnvcaller) {
        ch_svcaller_priority = Channel.value(["tiddit", "manta", "cnvnator"])
    } else {
        ch_svcaller_priority = Channel.value(["tiddit", "manta", "gcnvcaller", "cnvnator"])
    }

    // Read and store paths in the vep_plugin_files file
    if (params.vep_plugin_files) {
        ch_vep_extra_files_unsplit.splitCsv ( header:true )
            .map { row ->
                f = file(row.vep_files[0])
                if(f.isFile() || f.isDirectory()){
                    return [f]
                } else {
                    error("\nVep database file ${f} does not exist.")
                }
            }
            .collect()
            .set {ch_vep_extra_files}
    }

    // Read and store hgnc ids in a channel
    ch_vep_filters_scout_fmt
        .map { it -> CustomFunctions.parseHgncIds(it.text) }
        .mix (ch_vep_filters_std_fmt)
        .toList()
        .set {ch_hgnc_ids}

    // Input QC
    if (!params.skip_fastqc) {
        FASTQC (ch_samplesheet)
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
        ch_samplesheet,
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

    if (!params.skip_mt_subsample) {
        SUBSAMPLE_MT(
            ch_mapped.mt_bam_bai,
            params.mt_subsample_rd,
            params.mt_subsample_seed
        )
        ch_versions   = ch_versions.mix(SUBSAMPLE_MT.out.versions)
    }

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
        ch_intervals_y,
        Channel.value(params.ngsbits_samplegender_method)
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
        ch_genome_chrsizes,
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
        ch_foundin_header,
        Channel.value(params.sentieon_dnascope_pcr_indel_model)
    )
    ch_versions = ch_versions.mix(CALL_SNV.out.versions)

    //
    // VARIANT EVALUATION
    //
    if (params.run_rtgvcfeval) {
        VARIANT_EVALUATION (
            CALL_SNV.out.genome_vcf_tabix,
            ch_genome_fai,
            ch_rtg_truthvcfs,
            ch_sdf
        )
        ch_versions = ch_versions.mix(VARIANT_EVALUATION.out.versions)
    }

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
            ch_sv_dbs,
            ch_sv_bedpedbs,
            params.genome,
            params.vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta,
            ch_genome_dictionary,
            ch_vep_extra_files
        ).set { ch_sv_annotate }
        ch_versions = ch_versions.mix(ch_sv_annotate.versions)

        GENERATE_CLINICAL_SET_SV(
            ch_sv_annotate.vcf_ann,
            ch_hgnc_ids
        )
        ch_versions = ch_versions.mix(GENERATE_CLINICAL_SET_SV.out.versions)

        ANN_CSQ_PLI_SV (
            GENERATE_CLINICAL_SET_SV.out.vcf,
            ch_variant_consequences_sv
        )
        ch_versions = ch_versions.mix(ANN_CSQ_PLI_SV.out.versions)

        RANK_VARIANTS_SV (
            ANN_CSQ_PLI_SV.out.vcf_ann,
            ch_pedfile,
            ch_reduced_penetrance,
            ch_score_config_sv
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_SV.out.versions)

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
            ch_samples,
            ch_scatter_split_intervals,
            ch_vep_extra_files,
            ch_genome_chrsizes
        ).set { ch_snv_annotate }
        ch_versions = ch_versions.mix(ch_snv_annotate.versions)

        GENERATE_CLINICAL_SET_SNV(
            ch_snv_annotate.vcf_ann,
            ch_hgnc_ids
        )
        ch_versions = ch_versions.mix(GENERATE_CLINICAL_SET_SNV.out.versions)

        ANN_CSQ_PLI_SNV (
            GENERATE_CLINICAL_SET_SNV.out.vcf,
            ch_variant_consequences_snv
        )
        ch_versions = ch_versions.mix(ANN_CSQ_PLI_SNV.out.versions)

        RANK_VARIANTS_SNV (
            ANN_CSQ_PLI_SNV.out.vcf_ann,
            ch_pedfile,
            ch_reduced_penetrance,
            ch_score_config_snv
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_SNV.out.versions)

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
            ch_vep_extra_files
        ).set { ch_mt_annotate }
        ch_versions = ch_versions.mix(ch_mt_annotate.versions)

        GENERATE_CLINICAL_SET_MT(
            ch_mt_annotate.vcf_ann,
            ch_hgnc_ids
        )
        ch_versions = ch_versions.mix(GENERATE_CLINICAL_SET_MT.out.versions)

        ANN_CSQ_PLI_MT(
            GENERATE_CLINICAL_SET_MT.out.vcf,
            ch_variant_consequences_snv
        )
        ch_versions = ch_versions.mix(ANN_CSQ_PLI_MT.out.versions)

        RANK_VARIANTS_MT (
            ANN_CSQ_PLI_MT.out.vcf_ann,
            ch_pedfile,
            ch_reduced_penetrance,
            ch_score_config_mt
        )
        ch_versions = ch_versions.mix(RANK_VARIANTS_MT.out.versions)

    }

    // STEP 1.7: SMNCOPYNUMBERCALLER
    RENAME_BAM_FOR_SMNCALLER(ch_mapped.genome_marked_bam, "bam").output
        .collect{it}
        .toList()
        .set { ch_bam_list }

    RENAME_BAI_FOR_SMNCALLER(ch_mapped.genome_marked_bai, "bam.bai").output
        .collect{it}
        .toList()
        .set { ch_bai_list }

    ch_case_info
        .combine(ch_bam_list)
        .combine(ch_bai_list)
        .set { ch_bams_bais }

    SMNCOPYNUMBERCALLER (
        ch_bams_bais
    )
    ch_versions = ch_versions.mix(RENAME_BAM_FOR_SMNCALLER.out.versions)
    ch_versions = ch_versions.mix(RENAME_BAI_FOR_SMNCALLER.out.versions)
    ch_versions = ch_versions.mix(SMNCOPYNUMBERCALLER.out.versions)

    // ped correspondence, sex check, ancestry check
    if (!params.skip_peddy) {
        PEDDY (
            CALL_SNV.out.genome_vcf.join(CALL_SNV.out.genome_tabix, failOnMismatch:true, failOnDuplicate:true),
            ch_pedfile
        )
        ch_versions = ch_versions.mix(PEDDY.out.versions.first())
    }

    // Generate CGH files from sequencing data, turned off by default
    if ( !params.skip_vcf2cytosure && params.analysis_type != "wes" ) {
        GENERATE_CYTOSURE_FILES (
            ch_sv_annotate.vcf_ann,
            ch_sv_annotate.tbi,
            ch_mapped.genome_marked_bam,
            ch_sample_id_map,
            ch_vcf2cytosure_blacklist
        )
        ch_versions = ch_versions.mix(GENERATE_CYTOSURE_FILES.out.versions)
    }

    // GENS
    if ( !params.skip_gens && params.analysis_type != "wes" ) {
        GENS (
            ch_mapped.genome_bam_bai,
            CALL_SNV.out.genome_gvcf,
            ch_genome_fasta,
            ch_genome_fai,
            ch_gens_interval_list,
            ch_gens_pon_female,
            ch_gens_pon_male,
            ch_gens_gnomad_pos,
            ch_case_info,
            ch_genome_dictionary
        )
        ch_versions = ch_versions.mix(GENS.out.versions)
    }

    CALL_MOBILE_ELEMENTS(
        ch_mapped.genome_bam_bai,
        ch_genome_fasta,
        ch_genome_fai,
        ch_me_references,
        ch_case_info,
        params.genome
    )
    ch_versions = ch_versions.mix(CALL_MOBILE_ELEMENTS.out.versions)

    if (!params.skip_me_annotation) {
        ANNOTATE_MOBILE_ELEMENTS(
            CALL_MOBILE_ELEMENTS.out.vcf,
            ch_me_svdb_resources,
            ch_genome_fasta,
            ch_genome_dictionary,
            ch_vep_cache,
            ch_variant_consequences_sv,
            ch_hgnc_ids,
            params.genome,
            params.vep_cache_version,
            ch_vep_extra_files
        )
        ch_versions = ch_versions.mix(ANNOTATE_MOBILE_ELEMENTS.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.fromPath("$projectDir/docs/images/nf-core-raredisease_logo_light.png", checkIfExists: true)
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))
    if (!params.skip_fastqc) {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    }
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.multiple_metrics.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.hs_metrics.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.qualimap_results.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.global_dist.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.cov.map{it[1]}.collect().ifEmpty([]))

    if (!params.skip_peddy) {
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.ped.map{it[1]}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.csv.map{it[1]}.collect().ifEmpty([]))
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
