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
include { samplesheetToList       } from 'plugin/nf-schema'
include { CREATE_HGNCIDS_FILE     } from './modules/local/create_hgncids_file'
include { CREATE_PEDIGREE_FILE    } from './modules/local/create_pedigree_file'
include { channelFromPath         } from './subworkflows/local/utils_nfcore_raredisease_pipeline'
include { channelFromPathWithMeta } from './subworkflows/local/utils_nfcore_raredisease_pipeline'
include { channelFromSamplesheet  } from './subworkflows/local/utils_nfcore_raredisease_pipeline'
include { parseSkipList           } from './subworkflows/local/utils_nfcore_raredisease_pipeline'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_raredisease_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_raredisease_pipeline'
include { PREPARE_REFERENCES      } from './subworkflows/local/prepare_references'
include { RAREDISEASE             } from './workflows/raredisease'
include { SCATTER_GENOME          } from './subworkflows/local/scatter_genome'

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
    ch_alignments
    ch_case_info
    ch_reads
    ch_samples
    val_aligner
    val_analysis_type
    val_bwa
    val_bwamem2
    val_bwameme
    val_cadd_resources
    val_call_interval
    val_concatenate_snv_calls
    val_extract_alignments
    val_fai
    val_fasta
    val_gcnvcaller_model
    val_genome
    val_gens_gnomad_pos
    val_gens_interval_list
    val_gens_pon_female
    val_gens_pon_male
    val_gnomad_af
    val_gnomad_af_idx
    val_homoplasmy_af_threshold
    val_intervals_wgs
    val_intervals_y
    val_known_dbsnp
    val_known_dbsnp_tbi
    val_mbuffer_mem
    val_ml_model
    val_mobile_element_references
    val_mobile_element_svdb_annotations
    val_mt_aligner
    val_mt_fasta
    val_mt_subsample_approach
    val_mt_subsample_rd
    val_mt_subsample_seed
    val_multiqc_samples
    val_ngsbits_samplegender_method
    val_par_bed
    val_platform
    val_ploidy_model
    val_readcount_intervals
    val_reduced_penetrance
    val_rtg_truthvcfs
    val_run_mt_for_wes
    val_run_rtgvcfeval
    val_sambamba_regions
    val_sample_id_map
    val_samtools_sort_threads
    val_save_mapped_as_cram
    val_score_config_mt
    val_score_config_snv
    val_score_config_sv
    val_sdf
    val_sentieon_dnascope_pcr_indel_model
    val_sequence_dictionary
    val_skip_tools
    val_skip_subworkflows
    val_subdepth
    val_svdb_query_bedpedbs
    val_svdb_query_dbs
    val_target_bed
    val_variant_caller
    val_variant_catalog
    val_variant_consequences_snv
    val_variant_consequences_sv
    val_vcf2cytosure_blacklist
    val_vcfanno_extra_resources
    val_vcfanno_lua
    val_vcfanno_resources
    val_vcfanno_toml
    val_vep_cache_version
    val_vep_filters
    val_vep_filters_scout_fmt
    val_vep_plugin_files
    val_verifybamid_svd_bed
    val_verifybamid_svd_mu
    val_verifybamid_svd_ud
    val_vep_cache

    main:

    //
    // WORKFLOW: Run pipeline
    //

    ch_versions = channel.empty()

    PREPARE_REFERENCES (
        val_aligner,
        val_analysis_type,
        val_bwa,
        val_bwamem2,
        val_bwameme,
        val_fai,
        val_fasta,
        val_gnomad_af,
        val_gnomad_af_idx,
        val_known_dbsnp,
        val_known_dbsnp_tbi,
        val_mt_aligner,
        val_mt_fasta,
        val_run_mt_for_wes,
        val_run_rtgvcfeval,
        val_sdf,
        val_sequence_dictionary,
        val_target_bed,
        val_vcfanno_extra_resources,
        val_vep_cache
    )
    .set { ch_references }

    ch_bait_intervals           = ch_references.bait_intervals
    ch_dbsnp                    = ch_references.dbsnp
    ch_dbsnp_tbi                = ch_references.dbsnp_tbi
    ch_genome_bwaindex          = ch_references.genome_bwa_index
    ch_genome_bwamem2index      = ch_references.genome_bwamem2_index
    ch_genome_bwamemeindex      = ch_references.genome_bwameme_index
    ch_genome_chrsizes          = ch_references.genome_chrom_sizes
    ch_genome_fai               = ch_references.genome_fai
    ch_genome_fasta             = ch_references.genome_fasta
    ch_genome_dictionary        = ch_references.genome_dict
    ch_genome_hisat2index       = ch_references.genome_hisat2_index
    ch_gnomad_af                = ch_references.gnomad_af_idx
    ch_mt_bwaindex              = ch_references.mt_bwa_index
    ch_mt_bwamem2index          = ch_references.mt_bwamem2_index
    ch_mt_dictionary            = ch_references.mt_dict
    ch_mt_fai                   = ch_references.mt_fai
    ch_mt_fasta                 = ch_references.mt_fasta
    ch_mt_intervals             = ch_references.mt_intervals
    ch_mt_lastdb                = ch_references.mt_last_index
    ch_mtshift_backchain        = ch_references.mtshift_backchain
    ch_mtshift_bwaindex         = ch_references.mtshift_bwa_index
    ch_mtshift_bwamem2index     = ch_references.mtshift_bwamem2_index
    ch_mtshift_dictionary       = ch_references.mtshift_dict
    ch_mtshift_fai              = ch_references.mtshift_fai
    ch_mtshift_fasta            = ch_references.mtshift_fasta
    ch_mtshift_intervals        = ch_references.mtshift_intervals
    ch_sdf                      = ch_references.sdf
    ch_target_bed               = ch_references.target_bed
    ch_target_intervals         = ch_references.target_intervals
    ch_vcfanno_extra            = ch_references.vcfanno_extra
    ch_vep_cache                = ch_references.vep_resources
    ch_versions                 = ch_versions.mix(ch_references.versions)

    // Using channelFromPath helper (val_x ? channel.fromPath(val_x).collect() : channel.value([]))
    ch_cadd_resources           = channelFromPath(val_cadd_resources, true)
    ch_multiqc_samples          = channelFromPath(val_multiqc_samples, true)
    ch_reduced_penetrance       = channelFromPath(val_reduced_penetrance, true)
    ch_rtg_truthvcfs            = channelFromPath(val_rtg_truthvcfs, true)
    ch_score_config_mt          = channelFromPath(val_score_config_mt, true)
    ch_score_config_snv         = channelFromPath(val_score_config_snv, true)
    ch_score_config_sv          = channelFromPath(val_score_config_sv, true)
    ch_vcf2cytosure_blacklist   = channelFromPath(val_vcf2cytosure_blacklist, true)
    ch_vcfanno_lua              = channelFromPath(val_vcfanno_lua, true)
    ch_vcfanno_toml             = channelFromPath(val_vcfanno_toml, true)

    // Using channelFromPath helper (val_x ? channel.fromPath(val_x).collect() : channel.empty())
    ch_gens_gnomad_pos          = channelFromPath(val_gens_gnomad_pos)
    ch_gens_interval_list       = channelFromPath(val_gens_interval_list)
    ch_intervals_wgs            = channelFromPath(val_intervals_wgs)
    ch_intervals_y              = channelFromPath(val_intervals_y)
    ch_me_svdb_resources        = channelFromPath(val_mobile_element_svdb_annotations)
    ch_readcount_intervals      = channelFromPath(val_readcount_intervals)
    ch_svd_bed                  = channelFromPath(val_verifybamid_svd_bed)
    ch_svd_mu                   = channelFromPath(val_verifybamid_svd_mu)
    ch_svd_ud                   = channelFromPath(val_verifybamid_svd_ud)

    // Using channelFromPathWithMeta helper (with simpleName). If filepath is null, returns, [[:],[]]
    ch_call_interval            = channelFromPathWithMeta(val_call_interval, true)
    ch_ml_model                 = channelFromPathWithMeta(val_ml_model, true)
    ch_variant_catalog          = channelFromPathWithMeta(val_variant_catalog, true)
    ch_variant_consequences_snv = channelFromPathWithMeta(val_variant_consequences_snv, true)
    ch_variant_consequences_sv  = channelFromPathWithMeta(val_variant_consequences_sv, true)

    // Using channelFromPathWithMeta helper (with simpleName). If filepath is null, returns, empty channel
    ch_gens_pon_female          = channelFromPathWithMeta(val_gens_pon_female)
    ch_gens_pon_male            = channelFromPathWithMeta(val_gens_pon_male)
    ch_ploidy_model             = channelFromPathWithMeta(val_ploidy_model)

    // Using channelFromPathWithMeta helper. Returns either an empty channel or [[:],[]] or a channel with custom ID.
    ch_par_bed                  = channelFromPathWithMeta(val_par_bed, true, "par_bed")
    ch_sambamba_bed             = channelFromPathWithMeta(val_sambamba_regions, false, 'sambamba')
    ch_vep_filters_std_fmt      = channelFromPathWithMeta(val_vep_filters, false, 'standard')
    ch_vep_filters_scout_fmt    = channelFromPathWithMeta(val_vep_filters_scout_fmt, false, 'scout')

    // Using channelFromSamplesheet helper. Returns either an empty channel or validated channel.
    ch_me_references            = channelFromSamplesheet(val_mobile_element_references, "${projectDir}/assets/mobile_element_references_schema.json", false)
    ch_me_svdb_resources        = channelFromSamplesheet(val_mobile_element_svdb_annotations, "${projectDir}/assets/svdb_query_vcf_schema.json")
    ch_sample_id_map            = channelFromSamplesheet(val_sample_id_map, "${projectDir}/assets/sample_id_map.json")
    ch_svdb_bedpedbs            = channelFromSamplesheet(val_svdb_query_bedpedbs, "${projectDir}/assets/svdb_query_bedpe_schema.json")
    ch_svdb_dbs                 = channelFromSamplesheet(val_svdb_query_dbs, "${projectDir}/assets/svdb_query_vcf_schema.json")

    ch_cadd_header              = channel.fromPath("$projectDir/assets/cadd_to_vcf_header_-1.0-.txt", checkIfExists: true).collect()
    ch_foundin_header           = channel.fromPath("$projectDir/assets/foundin.hdr", checkIfExists: true).collect()
    ch_ngsbits_method           = channel.value(val_ngsbits_samplegender_method)
    ch_sentieon_pcr_indel_model = channel.value(val_sentieon_dnascope_pcr_indel_model)
    ch_subdepth                 = channel.value(val_subdepth)
    ch_vcfanno_resources        = val_vcfanno_resources ? channel.fromPath(val_vcfanno_resources).splitText().map{it -> it.trim()}.collect()
                                                        : channel.value([])
    ch_gcnvcaller_model         = val_gcnvcaller_model  ? channel.fromPath(val_gcnvcaller_model)
                                                            .splitCsv ( header:true )
                                                            .map { row ->
                                                                    return [[id:file(row.models).simpleName], row.models]
                                                                }
                                                            : channel.empty()

    //
    // Read and store paths in the vep_plugin_files file
    //
    ch_vep_extra_files = channel.empty()
    if (val_vep_plugin_files) {
        channel.fromPath(val_vep_plugin_files)
            .collect()
            .splitCsv ( header:true )
            .map { row ->
                def f = file(row.vep_files[0])
                if(f.isFile() || f.isDirectory()){
                    return [f]
                } else {
                    error("\nVep database file ${f} does not exist.")
                }
            }
            .collect()
            .set {ch_vep_extra_files}
    }

    //
    // Dump all HGNC ids in a file
    //
    ch_vep_filters_scout_fmt
        .mix (ch_vep_filters_std_fmt)
        .set {ch_vep_filters}

    CREATE_HGNCIDS_FILE(ch_vep_filters)
        .txt
        .set {ch_hgnc_ids}

    //
    // Generate pedigree file
    //
    ch_pedfile  = CREATE_PEDIGREE_FILE(ch_samples.toList()).ped

    // Tools
    skip_fastp                 = parseSkipList(val_skip_tools, 'fastp')
    skip_fastqc                = parseSkipList(val_skip_tools, 'fastqc')
    skip_gens                  = parseSkipList(val_skip_tools, 'gens')
    skip_germlinecnvcaller     = parseSkipList(val_skip_tools, 'germlinecnvcaller')
    skip_haplogrep3            = parseSkipList(val_skip_tools, 'haplogrep3')
    skip_ngsbits               = parseSkipList(val_skip_tools, 'ngsbits')
    skip_peddy                 = parseSkipList(val_skip_tools, 'peddy')
    skip_qualimap              = parseSkipList(val_skip_tools, 'qualimap')
    skip_smncopynumbercaller   = parseSkipList(val_skip_tools, 'smncopynumbercaller')
    skip_vcf2cytosure          = parseSkipList(val_skip_tools, 'vcf2cytosure')

    // Subworkflows
    skip_me_annotation         = parseSkipList(val_skip_subworkflows, 'me_annotation')
    skip_me_calling            = parseSkipList(val_skip_subworkflows, 'me_calling')
    skip_mt_annotation         = parseSkipList(val_skip_subworkflows, 'mt_annotation')
    skip_mt_subsample          = parseSkipList(val_skip_subworkflows, 'mt_subsample')
    skip_repeat_annotation     = parseSkipList(val_skip_subworkflows, 'repeat_annotation')
    skip_repeat_calling        = parseSkipList(val_skip_subworkflows, 'repeat_calling')
    skip_snv_annotation        = parseSkipList(val_skip_subworkflows, 'snv_annotation')
    skip_snv_calling           = parseSkipList(val_skip_subworkflows, 'snv_calling')
    skip_sv_annotation         = parseSkipList(val_skip_subworkflows, 'sv_annotation')
    skip_sv_calling            = parseSkipList(val_skip_subworkflows, 'sv_calling')
    skip_generate_clinical_set = parseSkipList(val_skip_subworkflows, 'generate_clinical_set')

    //
    // SV caller priority
    //
    if (skip_germlinecnvcaller) {
        if (val_analysis_type.equals("wgs")) {
            ch_svcaller_priority = channel.value(["tiddit", "manta", "cnvnator"])
        } else {
            ch_svcaller_priority = channel.value([])
        }
    } else {
        if (val_analysis_type.equals("wgs")) {
            ch_svcaller_priority = channel.value(["tiddit", "manta", "gcnvcaller", "cnvnator"])
        } else {
            ch_svcaller_priority = channel.value(["manta", "gcnvcaller"])
        }
    }

    //
    // Create chromosome bed and intervals for splitting and gathering operations
    //
    ch_scatter_split_intervals = channel.empty()
    if (!skip_snv_annotation) {
        SCATTER_GENOME (
            ch_genome_dictionary,
            ch_genome_fai,
            ch_genome_fasta
        ).split_intervals
        .set { ch_scatter_split_intervals }
    }

    RAREDISEASE (
        ch_alignments,
        ch_bait_intervals,
        ch_cadd_header,
        ch_cadd_resources,
        ch_call_interval,
        ch_case_info,
        ch_dbsnp,
        ch_dbsnp_tbi,
        ch_foundin_header,
        ch_gcnvcaller_model,
        ch_genome_bwaindex,
        ch_genome_bwamem2index,
        ch_genome_bwamemeindex,
        ch_genome_chrsizes,
        ch_genome_dictionary,
        ch_genome_fai,
        ch_genome_fasta,
        ch_genome_hisat2index,
        ch_gens_gnomad_pos,
        ch_gens_interval_list,
        ch_gens_pon_female,
        ch_gens_pon_male,
        ch_gnomad_af,
        ch_hgnc_ids,
        ch_intervals_wgs,
        ch_intervals_y,
        ch_me_references,
        ch_me_svdb_resources,
        ch_ml_model,
        ch_mt_bwaindex,
        ch_mt_bwamem2index,
        ch_mt_dictionary,
        ch_mt_fai,
        ch_mt_fasta,
        ch_mt_intervals,
        ch_mt_lastdb,
        ch_mtshift_backchain,
        ch_mtshift_bwaindex,
        ch_mtshift_bwamem2index,
        ch_mtshift_dictionary,
        ch_mtshift_fai,
        ch_mtshift_fasta,
        ch_mtshift_intervals,
        ch_multiqc_samples,
        ch_ngsbits_method,
        ch_par_bed,
        ch_pedfile,
        ch_ploidy_model,
        ch_readcount_intervals,
        ch_reads,
        ch_reduced_penetrance,
        ch_rtg_truthvcfs,
        ch_sambamba_bed,
        ch_sample_id_map,
        ch_samples,
        ch_scatter_split_intervals,
        ch_score_config_mt,
        ch_score_config_snv,
        ch_score_config_sv,
        ch_sdf,
        ch_sentieon_pcr_indel_model,
        ch_subdepth,
        ch_svcaller_priority,
        ch_svd_bed,
        ch_svd_mu,
        ch_svd_ud,
        ch_svdb_bedpedbs,
        ch_svdb_dbs,
        ch_target_bed,
        ch_target_intervals,
        ch_variant_catalog,
        ch_variant_consequences_snv,
        ch_variant_consequences_sv,
        ch_vcf2cytosure_blacklist,
        ch_vcfanno_extra,
        ch_vcfanno_lua,
        ch_vcfanno_resources,
        ch_vcfanno_toml,
        ch_vep_cache,
        ch_vep_extra_files,
        ch_versions,
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
        skip_fastp,
        skip_fastqc,
        skip_gens,
        skip_germlinecnvcaller,
        skip_haplogrep3,
        skip_ngsbits,
        skip_peddy,
        skip_qualimap,
        skip_smncopynumbercaller,
        skip_vcf2cytosure,
        val_aligner,
        val_analysis_type,
        val_cadd_resources,
        val_concatenate_snv_calls,
        val_extract_alignments,
        val_genome,
        val_homoplasmy_af_threshold,
        val_mbuffer_mem,
        val_mt_aligner,
        val_mt_subsample_approach,
        val_mt_subsample_rd,
        val_mt_subsample_seed,
        val_platform,
        val_run_mt_for_wes,
        val_run_rtgvcfeval,
        val_sample_id_map,
        val_samtools_sort_threads,
        val_save_mapped_as_cram,
        val_svdb_query_bedpedbs,
        val_svdb_query_dbs,
        val_target_bed,
        val_variant_caller,
        val_vep_cache_version
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
        PIPELINE_INITIALISATION.out.align,
        PIPELINE_INITIALISATION.out.case_info,
        PIPELINE_INITIALISATION.out.reads,
        PIPELINE_INITIALISATION.out.samples,
        params.aligner,
        params.analysis_type,
        params.bwa,
        params.bwamem2,
        params.bwameme,
        params.cadd_resources,
        params.call_interval,
        params.concatenate_snv_calls,
        params.extract_alignments,
        params.fai,
        params.fasta,
        params.gcnvcaller_model,
        params.genome,
        params.gens_gnomad_pos,
        params.gens_interval_list,
        params.gens_pon_female,
        params.gens_pon_male,
        params.gnomad_af,
        params.gnomad_af_idx,
        params.homoplasmy_af_threshold,
        params.intervals_wgs,
        params.intervals_y,
        params.known_dbsnp,
        params.known_dbsnp_tbi,
        params.mbuffer_mem,
        params.ml_model,
        params.mobile_element_references,
        params.mobile_element_svdb_annotations,
        params.mt_aligner,
        params.mt_fasta,
        params.mt_subsample_approach,
        params.mt_subsample_rd,
        params.mt_subsample_seed,
        params.multiqc_samples,
        params.ngsbits_samplegender_method,
        params.par_bed,
        params.platform,
        params.ploidy_model,
        params.readcount_intervals,
        params.reduced_penetrance,
        params.rtg_truthvcfs,
        params.run_mt_for_wes,
        params.run_rtgvcfeval,
        params.sambamba_regions,
        params.sample_id_map,
        params.samtools_sort_threads,
        params.save_mapped_as_cram,
        params.score_config_mt,
        params.score_config_snv,
        params.score_config_sv,
        params.sdf,
        params.sentieon_dnascope_pcr_indel_model,
        params.sequence_dictionary,
        params.skip_tools,
        params.skip_subworkflows,
        params.mitosalt_depth,
        params.svdb_query_bedpedbs,
        params.svdb_query_dbs,
        params.target_bed,
        params.variant_caller,
        params.variant_catalog,
        params.variant_consequences_snv,
        params.variant_consequences_sv,
        params.vcf2cytosure_blacklist,
        params.vcfanno_extra_resources,
        params.vcfanno_lua,
        params.vcfanno_resources,
        params.vcfanno_toml,
        params.vep_cache_version,
        params.vep_filters,
        params.vep_filters_scout_fmt,
        params.vep_plugin_files,
        params.verifybamid_svd_bed,
        params.verifybamid_svd_mu,
        params.verifybamid_svd_ud,
        params.vep_cache
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
