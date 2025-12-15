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

    main:

    //
    // WORKFLOW: Run pipeline
    //

    ch_versions                  = channel.empty()
    ch_genome_fasta              = channel.fromPath(params.fasta).map { it -> [[id:it.simpleName], it] }.collect()
    ch_gnomad_af_tab             = params.gnomad_af           ? channel.fromPath(params.gnomad_af).map{ it -> [[id:it.simpleName], it] }.collect()
                                                                : channel.value([[],[]])
    ch_dbsnp                     = params.known_dbsnp         ? channel.fromPath(params.known_dbsnp).map{ it -> [[id:it.simpleName], it] }.collect()
                                                                : channel.value([[],[]])
    ch_target_bed_unprocessed    = params.target_bed          ? channel.fromPath(params.target_bed).map{ it -> [[id:it.simpleName], it] }.collect()
                                                                : channel.value([[],[]])
    ch_vcfanno_extra_unprocessed = params.vcfanno_extra_resources ? channel.fromPath(params.vcfanno_extra_resources).map { it -> [[id:it.baseName], it] }.collect()
                                                                : channel.empty()
    ch_vep_cache_unprocessed     = params.vep_cache           ? channel.fromPath(params.vep_cache).map { it -> [[id:'vep_cache'], it] }.collect()
                                                                : channel.value([[],[]])

    PREPARE_REFERENCES (
        ch_genome_fasta,
        ch_gnomad_af_tab,
        ch_dbsnp,
        ch_target_bed_unprocessed,
        ch_vcfanno_extra_unprocessed,
        ch_vep_cache_unprocessed,
        params.aligner,
        params.analysis_type,
        params.bwa,
        params.bwamem2,
        params.bwameme,
        params.fai,
        params.mt_aligner,
        params.mt_fasta,
        params.run_mt_for_wes,
        params.sequence_dictionary
    )
    .set { ch_references }

    ch_bait_intervals           = ch_references.bait_intervals
    ch_cadd_header              = channel.fromPath("$projectDir/assets/cadd_to_vcf_header_-1.0-.txt", checkIfExists: true).collect()
    ch_cadd_resources           = params.cadd_resources                     ? channel.fromPath(params.cadd_resources).collect()
                                                                            : channel.value([])
    ch_call_interval            = params.call_interval                      ? channel.fromPath(params.call_interval).map {it -> [[id:it.simpleName], it]}.collect()
                                                                            : channel.value([[:],[]])
    ch_dbsnp_tbi                = params.known_dbsnp_tbi                    ? channel.fromPath(params.known_dbsnp_tbi).map {it -> [[id:it.simpleName], it]}.collect()
                                                                            : ch_references.known_dbsnp_tbi.ifEmpty([[],[]])
    ch_foundin_header           = channel.fromPath("$projectDir/assets/foundin.hdr", checkIfExists: true).collect()
    ch_gcnvcaller_model         = params.gcnvcaller_model                   ? channel.fromPath(params.gcnvcaller_model).splitCsv ( header:true )
                                                                            .map { row ->
                                                                                return [[id:file(row.models).simpleName], row.models]
                                                                            }
                                                                            : channel.empty()
    ch_genome_bwaindex          = params.bwa                                ? channel.fromPath(params.bwa).map {it -> [[id:it.simpleName], it]}.collect()
                                                                            : ch_references.genome_bwa_index
    ch_genome_bwamem2index      = params.bwamem2                            ? channel.fromPath(params.bwamem2).map {it -> [[id:it.simpleName], it]}.collect()
                                                                            : ch_references.genome_bwamem2_index
    ch_genome_bwamemeindex      = params.bwameme                            ? channel.fromPath(params.bwameme).map {it -> [[id:it.simpleName], it]}.collect()
                                                                            : ch_references.genome_bwameme_index
    ch_genome_chrsizes          = ch_references.genome_chrom_sizes
    ch_genome_fai               = ch_references.genome_fai
    ch_genome_dictionary        = ch_references.genome_dict
    ch_gens_gnomad_pos          = params.gens_gnomad_pos                    ? channel.fromPath(params.gens_gnomad_pos).collect()
                                                                            : channel.empty()
    ch_gens_interval_list       = params.gens_interval_list                 ? channel.fromPath(params.gens_interval_list).collect()
                                                                            : channel.empty()
    ch_gens_pon_female          = params.gens_pon_female                    ? channel.fromPath(params.gens_pon_female).map { it -> [ [id:it.simpleName], it ] }.collect()
                                                                            : channel.empty()
    ch_gens_pon_male            = params.gens_pon_male                      ? channel.fromPath(params.gens_pon_male).map { it -> [ [id:it.simpleName], it ] }.collect()
                                                                            : channel.empty()
    ch_gnomad_afidx             = params.gnomad_af_idx                      ? channel.fromPath(params.gnomad_af_idx).collect()
                                                                            : ch_references.gnomad_af_idx
    ch_gnomad_af                = params.gnomad_af                          ? ch_gnomad_af_tab.join(ch_gnomad_afidx).map {_meta, tab, idx -> [tab,idx]}.collect()
                                                                            : channel.empty()
    ch_intervals_wgs            = params.intervals_wgs                      ? channel.fromPath(params.intervals_wgs).collect()
                                                                            : channel.empty()
    ch_intervals_y              = params.intervals_y                        ? channel.fromPath(params.intervals_y).collect()
                                                                            : channel.empty()
    ch_me_references            = params.mobile_element_references          ? channel.fromList(samplesheetToList(params.mobile_element_references, "${projectDir}/assets/mobile_element_references_schema.json"))
                                                                            : channel.empty()
    ch_me_svdb_resources        = params.mobile_element_svdb_annotations    ? channel.fromPath(params.mobile_element_svdb_annotations)
                                                                            : channel.empty()
    ch_ml_model                 = params.variant_caller.equals("sentieon")  ? channel.fromPath(params.ml_model).map {it -> [[id:it.simpleName], it]}.collect()
                                                                            : channel.value([[:],[]])
    ch_mt_intervals             = ch_references.mt_intervals
    ch_mt_bwaindex              = ch_references.mt_bwa_index
    ch_mt_bwamem2index          = ch_references.mt_bwamem2_index
    ch_mt_dictionary            = ch_references.mt_dict
    ch_mt_fai                   = ch_references.mt_fai
    ch_mt_fasta                 = ch_references.mt_fasta
    ch_mtshift_backchain        = ch_references.mtshift_backchain
    ch_mtshift_bwaindex         = ch_references.mtshift_bwa_index
    ch_mtshift_bwamem2index     = ch_references.mtshift_bwamem2_index
    ch_mtshift_dictionary       = ch_references.mtshift_dict
    ch_mtshift_fai              = ch_references.mtshift_fai
    ch_mtshift_fasta            = ch_references.mtshift_fasta
    ch_mtshift_intervals        = ch_references.mtshift_intervals
    ch_par_bed                  = params.par_bed                            ? channel.fromPath(params.par_bed).map{ it -> [[id:'par_bed'], it] }.collect()
                                                                            : channel.value([[],[]])
    ch_ploidy_model             = params.ploidy_model                       ? channel.fromPath(params.ploidy_model).map{ it -> [[id:it.simpleName], it] }.collect()
                                                                            : channel.empty()
    ch_readcount_intervals      = params.readcount_intervals                ? channel.fromPath(params.readcount_intervals).collect()
                                                                            : channel.empty()
    ch_reduced_penetrance       = params.reduced_penetrance                 ? channel.fromPath(params.reduced_penetrance).collect()
                                                                            : channel.value([])
    ch_rtg_truthvcfs            = params.rtg_truthvcfs                      ? channel.fromPath(params.rtg_truthvcfs).collect()
                                                                            : channel.value([])
    ch_sambamba_bed             = params.sambamba_regions                   ? channel.fromPath(params.sambamba_regions).map{ it -> [[id:'sambamba'], it] }.collect()
                                                                            : channel.empty()
    ch_sample_id_map            = params.sample_id_map                      ? channel.fromList(samplesheetToList(params.sample_id_map, "${projectDir}/assets/sample_id_map.json"))
                                                                            : channel.empty()
    ch_score_config_mt          = params.score_config_mt                    ? channel.fromPath(params.score_config_mt).collect()
                                                                            : channel.value([])
    ch_score_config_snv         = params.score_config_snv                   ? channel.fromPath(params.score_config_snv).collect()
                                                                            : channel.value([])
    ch_score_config_sv          = params.score_config_sv                    ? channel.fromPath(params.score_config_sv).collect()
                                                                            : channel.value([])
    ch_sdf                      = params.sdf                                ? channel.fromPath(params.sdf).map{it -> [[id:it.simpleName],it]}.collect()
                                                                            : ch_references.sdf
    ch_svd_bed                  = params.verifybamid_svd_bed                ? channel.fromPath(params.verifybamid_svd_bed)
                                                                            : channel.empty()
    ch_svd_mu                   = params.verifybamid_svd_mu                 ? channel.fromPath(params.verifybamid_svd_mu)
                                                                            : channel.empty()
    ch_svd_ud                   = params.verifybamid_svd_ud                 ? channel.fromPath(params.verifybamid_svd_ud)
                                                                            : channel.empty()
    ch_svdb_bedpedbs            = params.svdb_query_bedpedbs                ? channel.fromList(samplesheetToList(params.svdb_query_bedpedbs, "assets/svdb_query_bedpe_schema.json")).collect()
                                                                            : channel.empty()
    ch_svdb_dbs                 = params.svdb_query_dbs                     ? channel.fromList(samplesheetToList(params.svdb_query_dbs, "assets/svdb_query_vcf_schema.json")).collect()
                                                                            : channel.empty()
    ch_target_bed               = ch_references.target_bed
    ch_target_intervals         = ch_references.target_intervals
    ch_variant_catalog          = params.variant_catalog                    ? channel.fromPath(params.variant_catalog).map { it -> [[id:it.simpleName],it]}.collect()
                                                                            : channel.value([[],[]])
    ch_variant_consequences_snv = params.variant_consequences_snv           ? channel.fromPath(params.variant_consequences_snv).map { it -> [[id:it.simpleName],it]}.collect()
                                                                            : channel.value([[],[]])
    ch_variant_consequences_sv  = params.variant_consequences_sv            ? channel.fromPath(params.variant_consequences_sv).map { it -> [[id:it.simpleName],it]}.collect()
                                                                            : channel.value([[],[]])
    ch_vcfanno_extra            = ch_references.vcfanno_extra
    ch_vcfanno_resources        = params.vcfanno_resources                  ? channel.fromPath(params.vcfanno_resources).splitText().map{it -> it.trim()}.collect()
                                                                            : channel.value([])
    ch_vcf2cytosure_blacklist   = params.vcf2cytosure_blacklist             ? channel.fromPath(params.vcf2cytosure_blacklist).collect()
                                                                            : channel.value([])
    ch_vcfanno_lua              = params.vcfanno_lua                        ? channel.fromPath(params.vcfanno_lua).collect()
                                                                            : channel.value([])
    ch_vcfanno_toml             = params.vcfanno_toml                       ? channel.fromPath(params.vcfanno_toml).collect()
                                                                            : channel.value([])
    ch_vep_cache                = ( params.vep_cache && params.vep_cache.endsWith("tar.gz") )  ? ch_references.vep_resources
                                                                            : ( params.vep_cache    ? channel.fromPath(params.vep_cache).collect() : channel.value([]) )
    ch_vep_filters_std_fmt      = params.vep_filters                        ? channel.fromPath(params.vep_filters).map { it -> [[id:'standard'],it]}.collect()
                                                                            : channel.empty()
    ch_vep_filters_scout_fmt    = params.vep_filters_scout_fmt              ? channel.fromPath(params.vep_filters_scout_fmt).map { it -> [[id:'scout'],it]}.collect()
                                                                            : channel.empty()
    ch_versions                 = ch_versions.mix(ch_references.versions)

    //
    // Read and store paths in the vep_plugin_files file
    //
    ch_vep_extra_files = channel.empty()
    if (params.vep_plugin_files) {
        channel.fromPath(params.vep_plugin_files)
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
    ch_versions = ch_versions.mix(CREATE_PEDIGREE_FILE.out.versions)

    // Tools
    skip_eklipse               = parseSkipList(params.skip_tools, 'eklipse')
    skip_fastp                 = parseSkipList(params.skip_tools, 'fastp')
    skip_fastqc                = parseSkipList(params.skip_tools, 'fastqc')
    skip_gens                  = parseSkipList(params.skip_tools, 'gens')
    skip_germlinecnvcaller     = parseSkipList(params.skip_tools, 'germlinecnvcaller')
    skip_haplogrep3            = parseSkipList(params.skip_tools, 'haplogrep3')
    skip_ngsbits               = parseSkipList(params.skip_tools, 'ngsbits')
    skip_peddy                 = parseSkipList(params.skip_tools, 'peddy')
    skip_qualimap              = parseSkipList(params.skip_tools, 'qualimap')
    skip_smncopynumbercaller   = parseSkipList(params.skip_tools, 'smncopynumbercaller')
    skip_vcf2cytosure          = parseSkipList(params.skip_tools, 'vcf2cytosure')

    // Subworkflows
    skip_me_annotation         = parseSkipList(params.skip_subworkflows, 'me_annotation')
    skip_me_calling            = parseSkipList(params.skip_subworkflows, 'me_calling')
    skip_mt_annotation         = parseSkipList(params.skip_subworkflows, 'mt_annotation')
    skip_mt_subsample          = parseSkipList(params.skip_subworkflows, 'mt_subsample')
    skip_repeat_annotation     = parseSkipList(params.skip_subworkflows, 'repeat_annotation')
    skip_repeat_calling        = parseSkipList(params.skip_subworkflows, 'repeat_calling')
    skip_snv_annotation        = parseSkipList(params.skip_subworkflows, 'snv_annotation')
    skip_snv_calling           = parseSkipList(params.skip_subworkflows, 'snv_calling')
    skip_sv_annotation         = parseSkipList(params.skip_subworkflows, 'sv_annotation')
    skip_sv_calling            = parseSkipList(params.skip_subworkflows, 'sv_calling')
    skip_generate_clinical_set = parseSkipList(params.skip_subworkflows, 'generate_clinical_set')

    //
    // SV caller priority
    //
    if (skip_germlinecnvcaller) {
        if (params.analysis_type.equals("wgs")) {
            ch_svcaller_priority = channel.value(["tiddit", "manta", "cnvnator"])
        } else {
            ch_svcaller_priority = channel.value([])
        }
    } else {
        if (params.analysis_type.equals("wgs")) {
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
        ch_mtshift_backchain,
        ch_mtshift_bwaindex,
        ch_mtshift_bwamem2index,
        ch_mtshift_dictionary,
        ch_mtshift_fai,
        ch_mtshift_fasta,
        ch_mtshift_intervals,
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
        params.analysis_type,
        params.svdb_query_bedpedbs,
        params.svdb_query_dbs,
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
        PIPELINE_INITIALISATION.out.samples
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
