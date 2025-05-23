/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap;samplesheetToList } from 'plugin/nf-schema'
include { paramsSummaryMultiqc               } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML             } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText             } from '../subworkflows/local/utils_nfcore_raredisease_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//
// MODULE: Installed directly from nf-core/modules
//

include { FASTQC                                            } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                           } from '../modules/nf-core/multiqc/main'
include { PEDDY                                             } from '../modules/nf-core/peddy/main'
include { SMNCOPYNUMBERCALLER                               } from '../modules/nf-core/smncopynumbercaller/main'
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_R1_FQ   } from '../modules/nf-core/spring/decompress/main'
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_R2_FQ   } from '../modules/nf-core/spring/decompress/main'
include { SPRING_DECOMPRESS as SPRING_DECOMPRESS_TO_FQ_PAIR } from '../modules/nf-core/spring/decompress/main'
include { STRANGER                                          } from '../modules/nf-core/stranger/main'

//
// MODULE: Local modules
//

include { RENAME_ALIGN_FILES as RENAME_BAM } from '../modules/local/rename_align_files'
include { RENAME_ALIGN_FILES as RENAME_BAI } from '../modules/local/rename_align_files'
include { CREATE_HGNCIDS_FILE              } from '../modules/local/create_hgncids_file'
include { CREATE_PEDIGREE_FILE             } from '../modules/local/create_pedigree_file'

//
// SUBWORKFLOWS
//

include { ALIGN                                              } from '../subworkflows/local/align'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_ME                 } from '../subworkflows/local/annotate_consequence_pli.nf'
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
include { GENERATE_CLINICAL_SET as GENERATE_CLINICAL_SET_ME  } from '../subworkflows/local/generate_clinical_set.nf'
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
    ch_reads
    ch_alignments
    ch_samples
    ch_case_info

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_mt_txt = Channel.empty()

    //
    // Initialize file channels for PREPARE_REFERENCES subworkflow
    //
    ch_genome_fasta              = Channel.fromPath(params.fasta).map { it -> [[id:it.simpleName], it] }.collect()
    ch_genome_fai                = params.fai                 ? Channel.fromPath(params.fai).map {it -> [[id:it.simpleName], it]}.collect()
                                                                : Channel.empty()
    ch_genome_dictionary         = params.sequence_dictionary ? Channel.fromPath(params.sequence_dictionary).map {it -> [[id:it.simpleName], it]}.collect()
                                                                : Channel.empty()
    ch_gnomad_af_tab             = params.gnomad_af           ? Channel.fromPath(params.gnomad_af).map{ it -> [[id:it.simpleName], it] }.collect()
                                                                : Channel.value([[],[]])
    ch_dbsnp                     = params.known_dbsnp         ? Channel.fromPath(params.known_dbsnp).map{ it -> [[id:it.simpleName], it] }.collect()
                                                                : Channel.value([[],[]])
    ch_mt_fasta                  = params.mt_fasta            ? Channel.fromPath(params.mt_fasta).map { it -> [[id:it.simpleName], it] }.collect()
                                                                : Channel.empty()
    ch_target_bed_unprocessed    = params.target_bed          ? Channel.fromPath(params.target_bed).map{ it -> [[id:it.simpleName], it] }.collect()
                                                                : Channel.value([[],[]])
    ch_vcfanno_extra_unprocessed = params.vcfanno_extra_resources ? Channel.fromPath(params.vcfanno_extra_resources).map { it -> [[id:it.baseName], it] }.collect()
                                                                : Channel.empty()
    ch_vep_cache_unprocessed     = params.vep_cache           ? Channel.fromPath(params.vep_cache).map { it -> [[id:'vep_cache'], it] }.collect()
                                                                : Channel.value([[],[]])

    //
    // Prepare references and indices.
    //
    PREPARE_REFERENCES (
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_dictionary,
        ch_mt_fasta,
        ch_gnomad_af_tab,
        ch_dbsnp,
        ch_target_bed_unprocessed,
        ch_vcfanno_extra_unprocessed,
        ch_vep_cache_unprocessed
    )
    .set { ch_references }

    //
    // Gather built indices or get them from the params
    //
    ch_bait_intervals           = ch_references.bait_intervals
    ch_cadd_header              = Channel.fromPath("$projectDir/assets/cadd_to_vcf_header_-1.0-.txt", checkIfExists: true).collect()
    ch_cadd_resources           = params.cadd_resources                     ? Channel.fromPath(params.cadd_resources).collect()
                                                                            : Channel.value([])
    ch_call_interval            = params.call_interval                      ? Channel.fromPath(params.call_interval).map {it -> [[id:it.simpleName], it]}.collect()
                                                                            : Channel.value([[:],[]])
    ch_dbsnp_tbi                = params.known_dbsnp_tbi                    ? Channel.fromPath(params.known_dbsnp_tbi).map {it -> [[id:it.simpleName], it]}.collect()
                                                                            : ch_references.known_dbsnp_tbi.ifEmpty([[],[]])
    ch_foundin_header           = Channel.fromPath("$projectDir/assets/foundin.hdr", checkIfExists: true).collect()
    ch_gcnvcaller_model         = params.gcnvcaller_model                   ? Channel.fromPath(params.gcnvcaller_model).splitCsv ( header:true )
                                                                            .map { row ->
                                                                                return [[id:file(row.models).simpleName], row.models]
                                                                            }
                                                                            : Channel.empty()
    ch_genome_bwaindex          = params.bwa                                ? Channel.fromPath(params.bwa).map {it -> [[id:it.simpleName], it]}.collect()
                                                                            : ch_references.genome_bwa_index
    ch_genome_bwamem2index      = params.bwamem2                            ? Channel.fromPath(params.bwamem2).map {it -> [[id:it.simpleName], it]}.collect()
                                                                            : ch_references.genome_bwamem2_index
    ch_genome_bwamemeindex      = params.bwameme                            ? Channel.fromPath(params.bwameme).map {it -> [[id:it.simpleName], it]}.collect()
                                                                            : ch_references.genome_bwameme_index
    ch_genome_chrsizes          = ch_references.genome_chrom_sizes
    ch_genome_fai               = ch_references.genome_fai
    ch_genome_dictionary        = ch_references.genome_dict
    ch_gens_gnomad_pos          = params.gens_gnomad_pos                    ? Channel.fromPath(params.gens_gnomad_pos).collect()
                                                                            : Channel.empty()
    ch_gens_interval_list       = params.gens_interval_list                 ? Channel.fromPath(params.gens_interval_list).collect()
                                                                            : Channel.empty()
    ch_gens_pon_female          = params.gens_pon_female                    ? Channel.fromPath(params.gens_pon_female).map { it -> [ [id:it.simpleName], it ] }.collect()
                                                                            : Channel.empty()
    ch_gens_pon_male            = params.gens_pon_male                      ? Channel.fromPath(params.gens_pon_male).map { it -> [ [id:it.simpleName], it ] }.collect()
                                                                            : Channel.empty()
    ch_gnomad_afidx             = params.gnomad_af_idx                      ? Channel.fromPath(params.gnomad_af_idx).collect()
                                                                            : ch_references.gnomad_af_idx
    ch_gnomad_af                = params.gnomad_af                          ? ch_gnomad_af_tab.join(ch_gnomad_afidx).map {meta, tab, idx -> [tab,idx]}.collect()
                                                                            : Channel.empty()
    ch_intervals_wgs            = params.intervals_wgs                      ? Channel.fromPath(params.intervals_wgs).collect()
                                                                            : Channel.empty()
    ch_intervals_y              = params.intervals_y                        ? Channel.fromPath(params.intervals_y).collect()
                                                                            : Channel.empty()
    ch_me_references            = params.mobile_element_references          ? Channel.fromList(samplesheetToList(params.mobile_element_references, "${projectDir}/assets/mobile_element_references_schema.json"))
                                                                            : Channel.empty()
    ch_me_svdb_resources        = params.mobile_element_svdb_annotations    ? Channel.fromPath(params.mobile_element_svdb_annotations)
                                                                            : Channel.empty()
    ch_ml_model                 = params.variant_caller.equals("sentieon")  ? Channel.fromPath(params.ml_model).map {it -> [[id:it.simpleName], it]}.collect()
                                                                            : Channel.value([[:],[]])
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
    ch_par_bed                  = params.par_bed                            ? Channel.fromPath(params.par_bed).map{ it -> [[id:'par_bed'], it] }.collect()
                                                                            : Channel.value([[],[]])
    ch_ploidy_model             = params.ploidy_model                       ? Channel.fromPath(params.ploidy_model).map{ it -> [[id:it.simpleName], it] }.collect()
                                                                            : Channel.empty()
    ch_readcount_intervals      = params.readcount_intervals                ? Channel.fromPath(params.readcount_intervals).collect()
                                                                            : Channel.empty()
    ch_reduced_penetrance       = params.reduced_penetrance                 ? Channel.fromPath(params.reduced_penetrance).collect()
                                                                            : Channel.value([])
    ch_rtg_truthvcfs            = params.rtg_truthvcfs                      ? Channel.fromPath(params.rtg_truthvcfs).collect()
                                                                            : Channel.value([])
    ch_sample_id_map            = params.sample_id_map                      ? Channel.fromList(samplesheetToList(params.sample_id_map, "${projectDir}/assets/sample_id_map.json"))
                                                                            : Channel.empty()
    ch_score_config_mt          = params.score_config_mt                    ? Channel.fromPath(params.score_config_mt).collect()
                                                                            : Channel.value([])
    ch_score_config_snv         = params.score_config_snv                   ? Channel.fromPath(params.score_config_snv).collect()
                                                                            : Channel.value([])
    ch_score_config_sv          = params.score_config_sv                    ? Channel.fromPath(params.score_config_sv).collect()
                                                                            : Channel.value([])
    ch_sdf                      = params.sdf                                ? Channel.fromPath(params.sdf).map{it -> [[id:it.simpleName],it]}.collect()
                                                                            : ch_references.sdf
    ch_sv_dbs                   = params.svdb_query_dbs                     ? Channel.fromPath(params.svdb_query_dbs)
                                                                            : Channel.empty()
    ch_sv_bedpedbs              = params.svdb_query_bedpedbs                ? Channel.fromPath(params.svdb_query_bedpedbs)
                                                                            : Channel.empty()
    ch_svd_bed                  = params.verifybamid_svd_bed                ? Channel.fromPath(params.verifybamid_svd_bed)
                                                                            : Channel.empty()
    ch_svd_mu                   = params.verifybamid_svd_mu                 ? Channel.fromPath(params.verifybamid_svd_mu)
                                                                            : Channel.empty()
    ch_svd_ud                   = params.verifybamid_svd_ud                 ? Channel.fromPath(params.verifybamid_svd_ud)
                                                                            : Channel.empty()
    ch_target_bed               = ch_references.target_bed
    ch_target_intervals         = ch_references.target_intervals
    ch_variant_catalog          = params.variant_catalog                    ? Channel.fromPath(params.variant_catalog).map { it -> [[id:it.simpleName],it]}.collect()
                                                                            : Channel.value([[],[]])
    ch_variant_consequences_snv = params.variant_consequences_snv           ? Channel.fromPath(params.variant_consequences_snv).collect()
                                                                            : Channel.value([])
    ch_variant_consequences_sv  = params.variant_consequences_sv            ? Channel.fromPath(params.variant_consequences_sv).collect()
                                                                            : Channel.value([])
    ch_vcfanno_extra            = ch_references.vcfanno_extra
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
    ch_vep_filters_std_fmt      = params.vep_filters                        ? Channel.fromPath(params.vep_filters).map { it -> [[id:'standard'],it]}.collect()
                                                                            : Channel.empty()
    ch_vep_filters_scout_fmt    = params.vep_filters_scout_fmt              ? Channel.fromPath(params.vep_filters_scout_fmt).map { it -> [[id:'scout'],it]}.collect()
                                                                            : Channel.empty()
    ch_versions                 = ch_versions.mix(ch_references.versions)

    //
    // SV caller priority
    //
    if (params.skip_tools && params.skip_tools.split(',').contains('germlinecnvcaller')) {
        if (params.analysis_type.equals("wgs")) {
            ch_svcaller_priority = Channel.value(["tiddit", "manta", "cnvnator"])
        } else {
            ch_svcaller_priority = Channel.value([])
        }
    } else {
        if (params.analysis_type.equals("wgs")) {
            ch_svcaller_priority = Channel.value(["tiddit", "manta", "gcnvcaller", "cnvnator"])
        } else {
            ch_svcaller_priority = Channel.value(["manta", "gcnvcaller"])
        }
    }

    //
    // Generate pedigree file
    //
    ch_pedfile   = CREATE_PEDIGREE_FILE(ch_samples.toList()).ped
    ch_versions = ch_versions.mix(CREATE_PEDIGREE_FILE.out.versions)

    //
    // Read and store paths in the vep_plugin_files file
    //
    ch_vep_extra_files = Channel.empty()
    if (params.vep_plugin_files) {
        ch_vep_extra_files_unsplit.splitCsv ( header:true )
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
    // Input QC (ch_reads will be empty if fastq input isn't provided so FASTQC won't run if input is not fastq)
    //

    ch_input_by_sample_type = ch_reads.branch{
        fastq_gz:           it[0].data_type == "fastq_gz"
        interleaved_spring: it[0].data_type == "interleaved_spring"
        separate_spring:    it[0].data_type == "separate_spring"
    }

    // Just one fastq.gz.spring-file with both R1 and R2
    ch_one_fastq_gz_pair_from_spring = SPRING_DECOMPRESS_TO_FQ_PAIR(ch_input_by_sample_type.interleaved_spring, false).fastq
    ch_versions                      = ch_versions.mix(SPRING_DECOMPRESS_TO_FQ_PAIR.out.versions.first())

    // Two fastq.gz.spring-files - one for R1 and one for R2
    ch_r1_fastq_gz_from_spring  = SPRING_DECOMPRESS_TO_R1_FQ(ch_input_by_sample_type.separate_spring.map{ meta, files -> [meta, files[0] ]}, true).fastq
    ch_r2_fastq_gz_from_spring  = SPRING_DECOMPRESS_TO_R2_FQ(ch_input_by_sample_type.separate_spring.map{ meta, files -> [meta, files[1] ]}, true).fastq
    ch_two_fastq_gz_from_spring = ch_r1_fastq_gz_from_spring.join(ch_r2_fastq_gz_from_spring).map{ meta, fastq_1, fastq_2 -> [meta, [fastq_1, fastq_2]]}
    ch_versions                 = ch_versions.mix(SPRING_DECOMPRESS_TO_R1_FQ.out.versions.first())
    ch_versions                 = ch_versions.mix(SPRING_DECOMPRESS_TO_R2_FQ.out.versions.first())

    ch_input_fastqs = ch_input_by_sample_type.fastq_gz.mix(ch_one_fastq_gz_pair_from_spring).mix(ch_two_fastq_gz_from_spring)

    //
    // Create chromosome bed and intervals for splitting and gathering operations
    //
    ch_scatter_split_intervals = Channel.empty()
    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('snv_annotation'))) {
        SCATTER_GENOME (
            ch_genome_dictionary,
            ch_genome_fai,
            ch_genome_fasta
        ).split_intervals
        .set { ch_scatter_split_intervals }
    }

    //
    // Input QC (ch_reads will be empty if fastq input isn't provided so FASTQC won't run if input is nott fastq)
    //

    if (!(params.skip_tools && params.skip_tools.split(',').contains('fastqc'))) {
        FASTQC (ch_input_fastqs)
        fastqc_report = FASTQC.out.zip
        ch_versions   = ch_versions.mix(FASTQC.out.versions.first())
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALIGN & FETCH STATS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ALIGN (
        ch_input_fastqs,
        ch_alignments,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_bwaindex,
        ch_genome_bwamem2index,
        ch_genome_bwamemeindex,
        ch_genome_dictionary,
        ch_mt_bwaindex,
        ch_mt_bwamem2index,
        ch_mt_dictionary,
        ch_mt_fai,
        ch_mt_fasta,
        ch_mtshift_bwaindex,
        ch_mtshift_bwamem2index,
        ch_mtshift_dictionary,
        ch_mtshift_fai,
        ch_mtshift_fasta,
        params.mbuffer_mem,
        params.platform,
        params.samtools_sort_threads
    )
    .set { ch_mapped }
    ch_versions   = ch_versions.mix(ALIGN.out.versions)

    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('mt_subsample')) && (params.analysis_type.equals("wgs") || params.run_mt_for_wes)) {
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
        ch_svd_bed,
        ch_svd_mu,
        ch_svd_ud,
        Channel.value(params.ngsbits_samplegender_method)
    )
    ch_versions = ch_versions.mix(QC_BAM.out.versions)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RENAME ALIGNMENT FILES FOR SMNCOPYNUMBERCALLER & REPEATCALLING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    if ( params.analysis_type.equals("wgs") && (!(params.skip_tools && params.skip_tools.split(',').contains('smncopynumbercaller')) || !(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('repeat_calling')))) {
        RENAME_BAM(ch_mapped.genome_marked_bam, "bam")
        RENAME_BAI(ch_mapped.genome_marked_bai, "bam.bai")
        ch_versions = ch_versions.mix(RENAME_BAM.out.versions)
        ch_versions = ch_versions.mix(RENAME_BAI.out.versions)
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL AND ANNOTATE REPEAT EXPANSIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('repeat_calling')) && params.analysis_type.equals("wgs") ) {
        CALL_REPEAT_EXPANSIONS (
            RENAME_BAM.out.output.join(RENAME_BAI.out.output, failOnMismatch:true, failOnDuplicate:true),
            ch_variant_catalog,
            ch_case_info,
            ch_genome_fasta,
            ch_genome_fai
        )
        ch_versions = ch_versions.mix(CALL_REPEAT_EXPANSIONS.out.versions)

        if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('repeat_annotation'))) {
            STRANGER (
                CALL_REPEAT_EXPANSIONS.out.vcf,
                ch_variant_catalog
            )
            ch_versions = ch_versions.mix(STRANGER.out.versions)
        }
    }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL AND ANNOTATE NUCLEAR AND MITOCHONDRIAL SNVs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('snv_calling'))) {
        CALL_SNV (
            ch_mapped.genome_bam_bai,
            ch_mapped.mt_bam_bai,
            ch_mapped.mtshift_bam_bai,
            ch_genome_chrsizes,
            ch_genome_fasta,
            ch_genome_fai,
            ch_genome_dictionary,
            ch_mt_intervals,
            ch_mt_dictionary,
            ch_mt_fai,
            ch_mt_fasta,
            ch_mtshift_dictionary,
            ch_mtshift_fai,
            ch_mtshift_fasta,
            ch_mtshift_intervals,
            ch_mtshift_backchain,
            ch_dbsnp,
            ch_dbsnp_tbi,
            ch_call_interval,
            ch_target_bed,
            ch_ml_model,
            ch_par_bed,
            ch_case_info,
            ch_foundin_header,
            Channel.value(params.sentieon_dnascope_pcr_indel_model)
        )
        ch_versions = ch_versions.mix(CALL_SNV.out.versions)
        ch_mt_txt = CALL_SNV.out.mt_txt

        //
        // ANNOTATE GENOME SNVs
        //
        if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('snv_annotation'))) {

            ANNOTATE_GENOME_SNVS (
                CALL_SNV.out.genome_vcf_tabix,
                params.analysis_type,
                ch_cadd_header,
                ch_cadd_resources,
                ch_vcfanno_extra,
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

            ch_snv_annotate.vcf_ann
                .multiMap { meta, vcf ->
                    clinical: [ meta + [ set: "clinical" ], vcf ]
                    research: [ meta + [ set: "research" ], vcf ]
                }
                .set { ch_clin_research_snv_vcf }

            ch_clinical_snv_vcf = Channel.empty()
            if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('generate_clinical_set'))) {
                GENERATE_CLINICAL_SET_SNV(
                    ch_clin_research_snv_vcf.clinical,
                    ch_hgnc_ids,
                    false
                )
                .vcf
                .set { ch_clinical_snv_vcf }
                ch_versions = ch_versions.mix(GENERATE_CLINICAL_SET_SNV.out.versions)
            }

            ch_ann_csq_snv_in = ch_clinical_snv_vcf.mix(ch_clin_research_snv_vcf.research)

            ANN_CSQ_PLI_SNV (
                ch_ann_csq_snv_in,
                ch_variant_consequences_snv
            )
            ch_versions = ch_versions.mix(ANN_CSQ_PLI_SNV.out.versions)

            ANN_CSQ_PLI_SNV.out.vcf_ann
                .filter { it ->
                    if (it[0].probands.size()==0) {
                        log.warn("Skipping nuclear SNV ranking since no affected samples are detected in the case")
                    }
                    it[0].probands.size()>0
                }
                .set {ch_ranksnv_nuclear_in}

            RANK_VARIANTS_SNV (
                ch_ranksnv_nuclear_in,
                ch_pedfile,
                ch_reduced_penetrance,
                ch_score_config_snv
            )
            ch_versions = ch_versions.mix(RANK_VARIANTS_SNV.out.versions)
        }

        //
        // ANNOTATE MT SNVs
        //
        if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('mt_annotation')) && (params.run_mt_for_wes || params.analysis_type.matches("wgs|mito"))) {

            ANNOTATE_MT_SNVS (
                CALL_SNV.out.mt_vcf,
                CALL_SNV.out.mt_tabix,
                ch_cadd_header,
                ch_cadd_resources,
                ch_genome_fasta,
                ch_vcfanno_extra,
                ch_vcfanno_lua,
                ch_vcfanno_resources,
                ch_vcfanno_toml,
                params.genome,
                params.vep_cache_version,
                ch_vep_cache,
                ch_vep_extra_files
            ).set { ch_mt_annotate }
            ch_versions = ch_versions.mix(ch_mt_annotate.versions)

            ch_mt_annotate.vcf_ann
                .multiMap { meta, vcf ->
                    clinical: [ meta + [ set: "clinical" ], vcf ]
                    research: [ meta + [ set: "research" ], vcf ]
                }
                .set { ch_clin_research_mt_vcf }

            ch_clinical_mtsnv_vcf = Channel.empty()
            if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('generate_clinical_set'))) {
                GENERATE_CLINICAL_SET_MT(
                    ch_clin_research_mt_vcf.clinical,
                    ch_hgnc_ids,
                    true
                )
                .vcf
                .set { ch_clinical_mtsnv_vcf }
                ch_versions = ch_versions.mix(GENERATE_CLINICAL_SET_MT.out.versions)
            }

            ch_ann_csq_mtsnv_in = ch_clinical_mtsnv_vcf.mix(ch_clin_research_mt_vcf.research)

            ANN_CSQ_PLI_MT(
                ch_ann_csq_mtsnv_in,
                ch_variant_consequences_snv
            )
            ch_versions = ch_versions.mix(ANN_CSQ_PLI_MT.out.versions)

            ANN_CSQ_PLI_MT.out.vcf_ann
                .filter { it ->
                    if (it[0].probands.size()==0) {
                        log.warn("Skipping mitochondrial SNV ranking since no affected samples are detected in the case")
                    }
                    it[0].probands.size()>0
                }
                .set {ch_ranksnv_mt_in}

            RANK_VARIANTS_MT (
                ch_ranksnv_mt_in,
                ch_pedfile,
                ch_reduced_penetrance,
                ch_score_config_mt
            )
            ch_versions = ch_versions.mix(RANK_VARIANTS_MT.out.versions)
        }
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL AND ANNOTATE NUCLEAR SVs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('sv_calling'))) {
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
        if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('sv_annotation'))) {
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

            ch_sv_annotate.vcf_ann
                .multiMap { meta, vcf ->
                    clinical: [ meta + [ set: "clinical" ], vcf ]
                    research: [ meta + [ set: "research" ], vcf ]
                }
                .set { ch_clin_research_sv_vcf }

            ch_clinical_sv_vcf = Channel.empty()
            if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('generate_clinical_set'))) {
                GENERATE_CLINICAL_SET_SV(
                    ch_clin_research_sv_vcf.clinical,
                    ch_hgnc_ids,
                    false
                )
                .vcf
                .set { ch_clinical_sv_vcf }
                ch_versions = ch_versions.mix(GENERATE_CLINICAL_SET_SV.out.versions)
            }

            ch_ann_csq_sv_in = ch_clinical_sv_vcf.mix(ch_clin_research_sv_vcf.research)

            ANN_CSQ_PLI_SV (
                ch_ann_csq_sv_in,
                ch_variant_consequences_sv
            )
            ch_versions = ch_versions.mix(ANN_CSQ_PLI_SV.out.versions)

            ANN_CSQ_PLI_SV.out.vcf_ann
                .filter { it ->
                    if (it[0].probands.size()==0) {
                        log.warn("Skipping SV ranking since no affected samples are detected in the case")
                    }
                    it[0].probands.size()>0
                }
                .set {ch_ranksnv_sv_in}

            RANK_VARIANTS_SV (
                ch_ranksnv_sv_in,
                ch_pedfile,
                ch_reduced_penetrance,
                ch_score_config_sv
            )
            ch_versions = ch_versions.mix(RANK_VARIANTS_SV.out.versions)
        }
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL AND ANNOTATE MOBILE ELEMENTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('me_calling')) && params.analysis_type.equals("wgs")) {
        CALL_MOBILE_ELEMENTS(
            ch_mapped.genome_bam_bai,
            ch_genome_fasta,
            ch_genome_fai,
            ch_me_references,
            ch_case_info,
            params.genome
        )
        ch_versions = ch_versions.mix(CALL_MOBILE_ELEMENTS.out.versions)

        if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('me_annotation'))) {
            ANNOTATE_MOBILE_ELEMENTS(
                CALL_MOBILE_ELEMENTS.out.vcf,
                ch_me_svdb_resources,
                ch_genome_fasta,
                ch_genome_dictionary,
                ch_vep_cache,
                params.genome,
                params.vep_cache_version,
                ch_vep_extra_files
            ).set { ch_me_annotate }
            ch_versions = ch_versions.mix(ANNOTATE_MOBILE_ELEMENTS.out.versions)

            ch_me_annotate.vcf_ann
                .multiMap { meta, vcf ->
                    clinical: [ meta + [ set: "clinical" ], vcf ]
                    research: [ meta + [ set: "research" ], vcf ]
                }
                .set { ch_clin_research_me_vcf }

            ch_clinical_me_vcf = Channel.empty()
            if (!(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('generate_clinical_set'))) {
                GENERATE_CLINICAL_SET_ME(
                    ch_clin_research_me_vcf.clinical,
                    ch_hgnc_ids,
                    false
                )
                .vcf
                .set { ch_clinical_me_vcf }
                ch_versions = ch_versions.mix( GENERATE_CLINICAL_SET_ME.out.versions )
            }

            ch_ann_csq_me_in = ch_clinical_me_vcf.mix(ch_clin_research_me_vcf.research)

            ANN_CSQ_PLI_ME(
                ch_ann_csq_me_in,
                ch_variant_consequences_sv
            )
            ch_versions = ch_versions.mix( ANN_CSQ_PLI_ME.out.versions )

        }
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SMNCOPYNUMBERCALLER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    if ( params.analysis_type.equals("wgs") && !(params.skip_tools && params.skip_tools.split(',').contains('smncopynumbercaller')) ) {

        RENAME_BAM.out.output
            .collect{it[1]}
            .toList()
            .set { ch_bam_list }

        RENAME_BAI.out.output
            .collect{it[1]}
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
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PEDDY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    if (!(params.skip_tools && params.skip_tools.split(',').contains('peddy'))) {
        PEDDY (
            CALL_SNV.out.genome_vcf.join(CALL_SNV.out.genome_tabix, failOnMismatch:true, failOnDuplicate:true),
            ch_pedfile
        )
        ch_versions = ch_versions.mix(PEDDY.out.versions.first())
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Generate CGH files from sequencing data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    if ( !(params.skip_tools && params.skip_tools.split(',').contains('vcf2cytosure')) && params.analysis_type.equals("wgs") && !(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('sv_calling')) && !(params.skip_subworkflows && params.skip_subworkflows.split(',').contains('sv_annotation'))) {
        GENERATE_CYTOSURE_FILES (
            ch_sv_annotate.vcf_ann,
            ch_sv_annotate.tbi,
            ch_mapped.genome_marked_bam,
            ch_sample_id_map,
            ch_vcf2cytosure_blacklist
        )
        ch_versions = ch_versions.mix(GENERATE_CYTOSURE_FILES.out.versions)
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    if (!(params.skip_tools && params.skip_tools.split(',').contains('gens')) && params.analysis_type.equals("wgs")) {
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VARIANT EVALUATION WITH RTGTOOLS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    if (params.run_rtgvcfeval) {
        VARIANT_EVALUATION (
            CALL_SNV.out.genome_vcf_tabix,
            ch_genome_fai,
            ch_rtg_truthvcfs,
            ch_sdf
        )
        ch_versions = ch_versions.mix(VARIANT_EVALUATION.out.versions)
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT SOFTWARE VERSIONS & MultiQC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'raredisease_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.fromPath("$projectDir/docs/images/nf-core-raredisease_logo_light.png", checkIfExists: true)


    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    ch_multiqc_files = ch_multiqc_files.mix(fastqc_report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_mt_txt.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN.out.fastp_json.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN.out.markdup_metrics.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.sex_check.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.multiple_metrics.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.hs_metrics.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.qualimap_results.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.global_dist.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.cov.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.self_sm.map{it[1]}.collect().ifEmpty([]))

    if (!(params.skip_tools && params.skip_tools.split(',').contains('peddy'))) {
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.ped.map{it[1]}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.csv.map{it[1]}.collect().ifEmpty([]))
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
