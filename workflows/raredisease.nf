/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_raredisease_pipeline'

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

//
// SUBWORKFLOWS
//

include { ALIGN                                                       } from '../subworkflows/local/align'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_ME                          } from '../subworkflows/local/annotate_consequence_pli.nf'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_MT                          } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SNV                         } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_SV                          } from '../subworkflows/local/annotate_consequence_pli'
include { ANNOTATE_GENOME_SNVS                                        } from '../subworkflows/local/annotate_genome_snvs'
include { ANNOTATE_MOBILE_ELEMENTS                                    } from '../subworkflows/local/annotate_mobile_elements'
include { ANNOTATE_MT_SNVS                                            } from '../subworkflows/local/annotate_mt_snvs'
include { ANNOTATE_STRUCTURAL_VARIANTS                                } from '../subworkflows/local/annotate_structural_variants'
include { CALL_MOBILE_ELEMENTS                                        } from '../subworkflows/local/call_mobile_elements'
include { CALL_REPEAT_EXPANSIONS                                      } from '../subworkflows/local/call_repeat_expansions'
include { CALL_SNV                                                    } from '../subworkflows/local/call_snv'
include { CALL_STRUCTURAL_VARIANTS                                    } from '../subworkflows/local/call_structural_variants'
include { CALL_SV_MT                                                  } from '../subworkflows/local/call_sv_MT'
include { GENERATE_CYTOSURE_FILES                                     } from '../subworkflows/local/generate_cytosure_files'
include { GENS                                                        } from '../subworkflows/local/gens'
include { PREPARE_REFERENCES                                          } from '../subworkflows/local/prepare_references'
include { QC_BAM                                                      } from '../subworkflows/local/qc_bam'
include { RANK_VARIANTS as RANK_VARIANTS_MT                           } from '../subworkflows/local/rank_variants'
include { RANK_VARIANTS as RANK_VARIANTS_SNV                          } from '../subworkflows/local/rank_variants'
include { RANK_VARIANTS as RANK_VARIANTS_SV                           } from '../subworkflows/local/rank_variants'
include { SUBSAMPLE_MT_FRAC                                           } from '../subworkflows/local/subsample_mt_frac'
include { SUBSAMPLE_MT_READS                                          } from '../subworkflows/local/subsample_mt_reads'
include { VARIANT_EVALUATION                                          } from '../subworkflows/local/variant_evaluation'
include { VCF_FILTER_BCFTOOLS_ENSEMBLVEP as GENERATE_CLINICAL_SET_ME  } from '../subworkflows/nf-core/vcf_filter_bcftools_ensemblvep/main'
include { VCF_FILTER_BCFTOOLS_ENSEMBLVEP as GENERATE_CLINICAL_SET_MT  } from '../subworkflows/nf-core/vcf_filter_bcftools_ensemblvep/main'
include { VCF_FILTER_BCFTOOLS_ENSEMBLVEP as GENERATE_CLINICAL_SET_SNV } from '../subworkflows/nf-core/vcf_filter_bcftools_ensemblvep/main'
include { VCF_FILTER_BCFTOOLS_ENSEMBLVEP as GENERATE_CLINICAL_SET_SV  } from '../subworkflows/nf-core/vcf_filter_bcftools_ensemblvep/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RAREDISEASE {

    take:
    ch_alignments
    ch_bait_intervals
    ch_cadd_header
    ch_cadd_resources
    ch_call_interval
    ch_case_info
    ch_dbsnp
    ch_dbsnp_tbi
    ch_foundin_header
    ch_gcnvcaller_model
    ch_genome_bwaindex
    ch_genome_bwamem2index
    ch_genome_bwamemeindex
    ch_genome_chrsizes
    ch_genome_dictionary
    ch_genome_fai
    ch_genome_fasta
    ch_genome_hisat2index
    ch_gens_gnomad_pos
    ch_gens_interval_list
    ch_gens_pon_female
    ch_gens_pon_male
    ch_gnomad_af
    ch_hgnc_ids
    ch_intervals_wgs
    ch_intervals_y
    ch_me_references
    ch_me_svdb_resources
    ch_ml_model
    ch_mt_bwaindex
    ch_mt_bwamem2index
    ch_mt_dictionary
    ch_mt_fai
    ch_mt_fasta
    ch_mt_intervals
    ch_mt_lastdb
    ch_mtshift_backchain
    ch_mtshift_bwaindex
    ch_mtshift_bwamem2index
    ch_mtshift_dictionary
    ch_mtshift_fai
    ch_mtshift_fasta
    ch_mtshift_intervals
    ch_multiqc_samples
    ch_ngsbits_method
    ch_par_bed
    ch_pedfile
    ch_ploidy_model
    ch_readcount_intervals
    ch_reads
    ch_reduced_penetrance
    ch_rtg_truthvcfs
    ch_sambamba_bed
    ch_sample_id_map
    ch_samples
    ch_scatter_split_intervals
    ch_score_config_mt
    ch_score_config_snv
    ch_score_config_sv
    ch_sdf
    ch_sentieon_pcr_indel_model
    ch_subdepth
    ch_svcaller_priority
    ch_svd_bed
    ch_svd_mu
    ch_svd_ud
    ch_svdb_bedpedbs
    ch_svdb_dbs
    ch_target_bed
    ch_target_intervals
    ch_variant_catalog
    ch_variant_consequences_snv
    ch_variant_consequences_sv
    ch_vcf2cytosure_blacklist
    ch_vcfanno_extra
    ch_vcfanno_lua
    ch_vcfanno_resources
    ch_vcfanno_toml
    ch_vep_cache
    ch_vep_extra_files
    ch_versions
    skip_me_calling
    skip_me_annotation
    skip_mt_annotation
    skip_mt_subsample
    skip_repeat_annotation
    skip_repeat_calling
    skip_snv_annotation
    skip_snv_calling
    skip_sv_annotation
    skip_sv_calling
    skip_generate_clinical_set
    skip_fastp
    skip_fastqc
    skip_gens
    skip_germlinecnvcaller
    skip_haplogrep3
    skip_ngsbits
    skip_peddy
    skip_qualimap
    skip_smncopynumbercaller
    skip_vcf2cytosure
    val_aligner
    val_analysis_type
    val_cadd_resources
    val_concatenate_snv_calls
    val_extract_alignments
    val_genome
    val_homoplasmy_af_threshold
    val_mbuffer_mem
    val_mt_aligner
    val_mt_subsample_approach
    val_mt_subsample_rd
    val_mt_subsample_seed
    val_platform
    val_run_mt_for_wes
    val_run_rtgvcfeval
    val_sample_id_map
    val_samtools_sort_threads
    val_save_mapped_as_cram
    val_svdb_query_bedpedbs
    val_svdb_query_dbs
    val_target_bed
    val_variant_caller
    val_vep_cache_version

    main:

    ch_multiqc_files = channel.empty()

    //
    // Input QC (ch_reads will be empty if fastq input isn't provided so FASTQC won't run if input is not fastq)
    //

    ch_input_by_sample_type = ch_reads.branch{ meta, _reads ->
        fastq_gz:           meta.data_type == "fastq_gz"
        interleaved_spring: meta.data_type == "interleaved_spring"
        separate_spring:    meta.data_type == "separate_spring"
    }

    // Just one fastq.gz.spring-file with both R1 and R2
    ch_one_fastq_gz_pair_from_spring = SPRING_DECOMPRESS_TO_FQ_PAIR(ch_input_by_sample_type.interleaved_spring, false).fastq

    // Two fastq.gz.spring-files - one for R1 and one for R2
    ch_r1_fastq_gz_from_spring  = SPRING_DECOMPRESS_TO_R1_FQ(ch_input_by_sample_type.separate_spring.map{ meta, files -> [meta, files[0] ]}, true).fastq
    ch_r2_fastq_gz_from_spring  = SPRING_DECOMPRESS_TO_R2_FQ(ch_input_by_sample_type.separate_spring.map{ meta, files -> [meta, files[1] ]}, true).fastq
    ch_two_fastq_gz_from_spring = ch_r1_fastq_gz_from_spring.join(ch_r2_fastq_gz_from_spring).map{ meta, fastq_1, fastq_2 -> [meta, [fastq_1, fastq_2]]}

    ch_input_fastqs = ch_input_by_sample_type.fastq_gz.mix(ch_one_fastq_gz_pair_from_spring).mix(ch_two_fastq_gz_from_spring)

    //
    // Input QC (ch_reads will be empty if fastq input isn't provided so FASTQC won't run if input is nott fastq)
    //

    if (!skip_fastqc) {
        FASTQC (ch_input_fastqs)
        fastqc_report = FASTQC.out.zip
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALIGN & FETCH STATS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ALIGN (
        ch_alignments,
        ch_genome_bwaindex,
        ch_genome_bwamem2index,
        ch_genome_bwamemeindex,
        ch_genome_dictionary,
        ch_genome_fai,
        ch_genome_fasta,
        ch_input_fastqs,
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
        skip_fastp,
        val_aligner,
        val_analysis_type,
        val_extract_alignments,
        val_mbuffer_mem,
        val_mt_aligner,
        val_platform,
        val_run_mt_for_wes,
        val_samtools_sort_threads,
        val_save_mapped_as_cram
    )
    .set { ch_mapped }

    if (!(skip_mt_subsample) && (val_analysis_type.equals("wgs") || val_run_mt_for_wes)) {
        if (val_mt_subsample_approach.equals("fraction")) {
            SUBSAMPLE_MT_FRAC(
                ch_mapped.mt_bam_bai,
                val_mt_subsample_rd,
                val_mt_subsample_seed
            )
            ch_versions   = ch_versions.mix(SUBSAMPLE_MT_FRAC.out.versions)
        } else {
            SUBSAMPLE_MT_READS(
                ch_mapped.mt_bam_bai,
            )
        }
    }

    //
    // BAM QUALITY CHECK
    //
    QC_BAM (
        ch_mapped.genome_marked_bam,
        ch_mapped.genome_marked_bam_bai,
        ch_bait_intervals,
        ch_genome_chrsizes,
        ch_genome_fai,
        ch_genome_fasta,
        ch_intervals_wgs,
        ch_intervals_y,
        ch_ngsbits_method,
        ch_svd_bed,
        ch_svd_mu,
        ch_svd_ud,
        ch_sambamba_bed,
        ch_target_intervals,
        val_analysis_type,
        val_aligner,
        val_target_bed,
        skip_ngsbits,
        skip_qualimap
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RENAME ALIGNMENT FILES FOR SMNCOPYNUMBERCALLER & REPEATCALLING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    if ( val_analysis_type.equals("wgs") && (!skip_smncopynumbercaller || !skip_repeat_calling)) {
        RENAME_BAM(ch_mapped.genome_marked_bam, "bam")
        RENAME_BAI(ch_mapped.genome_marked_bai, "bam.bai")
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL AND ANNOTATE REPEAT EXPANSIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    if (!skip_repeat_calling && val_analysis_type.equals("wgs") ) {
        CALL_REPEAT_EXPANSIONS (
            RENAME_BAM.out.output.join(RENAME_BAI.out.output, failOnMismatch:true, failOnDuplicate:true),
            ch_variant_catalog,
            ch_case_info,
            ch_genome_fasta,
            ch_genome_fai
        )

        if (!skip_repeat_annotation) {
            STRANGER (
                CALL_REPEAT_EXPANSIONS.out.vcf,
                ch_variant_catalog
            )
        }
    }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL AND ANNOTATE NUCLEAR AND MITOCHONDRIAL SNVs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    if (!skip_snv_calling) {
        CALL_SNV (
            ch_call_interval,
            ch_case_info,
            ch_dbsnp,
            ch_dbsnp_tbi,
            ch_foundin_header,
            ch_mapped.genome_marked_bam_bai,
            ch_genome_chrsizes,
            ch_genome_fasta,
            ch_genome_fai,
            ch_ml_model,
            ch_mapped.mt_bam_bai_gatksubwf,
            ch_mt_dictionary,
            ch_mt_fai,
            ch_mt_fasta,
            ch_mt_intervals,
            ch_mapped.mtshift_bam_bai_gatksubwf,
            ch_mtshift_dictionary,
            ch_mtshift_fai,
            ch_mtshift_fasta,
            ch_mtshift_intervals,
            ch_mtshift_backchain,
            ch_par_bed,
            ch_sentieon_pcr_indel_model,
            ch_target_bed,
            val_analysis_type,
            val_concatenate_snv_calls,
            val_run_mt_for_wes,
            val_variant_caller
        )
        ch_versions = ch_versions.mix(CALL_SNV.out.versions)

        //
        // ANNOTATE GENOME SNVs
        //
        if (!skip_snv_annotation) {

            ANNOTATE_GENOME_SNVS (
                ch_cadd_header,
                ch_cadd_resources,
                ch_genome_chrsizes,
                ch_genome_fai,
                ch_genome_fasta,
                ch_gnomad_af,
                ch_samples,
                ch_scatter_split_intervals,
                CALL_SNV.out.genome_vcf_tabix,
                ch_vcfanno_extra,
                ch_vcfanno_lua,
                ch_vcfanno_resources,
                ch_vcfanno_toml,
                ch_vep_cache,
                ch_vep_extra_files,
                val_cadd_resources,
                val_genome,
                val_vep_cache_version
            ).set { ch_snv_annotate }
            ch_versions = ch_versions.mix(ch_snv_annotate.versions)

            ch_snv_annotate.vcf_ann
                .multiMap { meta, vcf ->
                    clinical: [ meta + [ set: "clinical" ], vcf ]
                    research: [ meta + [ set: "research" ], vcf ]
                }
                .set { ch_clin_research_snv_vcf }

            ch_clinical_snv_vcf = channel.empty()
            if (!skip_generate_clinical_set) {
                GENERATE_CLINICAL_SET_SNV(
                    ch_clin_research_snv_vcf.clinical,
                    ch_hgnc_ids,
                    false,
                    true
                )
                GENERATE_CLINICAL_SET_SNV.out.vcf
                .set { ch_clinical_snv_vcf }
            }

            ch_ann_csq_snv_in = ch_clinical_snv_vcf.mix(ch_clin_research_snv_vcf.research)

            ANN_CSQ_PLI_SNV (
                ch_variant_consequences_snv,
                ch_ann_csq_snv_in,
                false
            )

            ANN_CSQ_PLI_SNV.out.vcf_ann
                .filter { meta, _vcf ->
                    if (meta.probands.size()==0) {
                        log.warn("Skipping nuclear SNV ranking since no affected samples are detected in the case")
                    }
                    meta.probands.size()>0
                }
                .set {ch_ranksnv_nuclear_in}

            RANK_VARIANTS_SNV (
                ch_pedfile,
                ch_reduced_penetrance,
                ch_score_config_snv,
                ch_ranksnv_nuclear_in,
                false
            )
        }

        //
        // ANNOTATE MT SNVs
        //
        if (!(skip_mt_annotation) && (val_run_mt_for_wes || val_analysis_type.matches("wgs|mito"))) {

            ANNOTATE_MT_SNVS (
                ch_cadd_header,
                ch_cadd_resources,
                ch_genome_fasta,
                ch_genome_fai,
                CALL_SNV.out.mt_vcf,
                ch_vcfanno_extra,
                ch_vcfanno_lua,
                ch_vcfanno_resources,
                ch_vcfanno_toml,
                ch_vep_cache,
                ch_vep_extra_files,
                skip_haplogrep3,
                val_cadd_resources,
                val_genome,
                val_homoplasmy_af_threshold,
                val_vep_cache_version
            ).set { ch_mt_annotate }
            ch_versions = ch_versions.mix(ch_mt_annotate.versions)

            ch_mt_annotate.vcf_ann
                .multiMap { meta, vcf ->
                    clinical: [ meta + [ set: "clinical" ], vcf ]
                    research: [ meta + [ set: "research" ], vcf ]
                }
                .set { ch_clin_research_mt_vcf }

            ch_clinical_mtsnv_vcf = channel.empty()
            if (!skip_generate_clinical_set) {
                GENERATE_CLINICAL_SET_MT(
                    ch_clin_research_mt_vcf.clinical,
                    ch_hgnc_ids,
                    true,
                    false
                )
                GENERATE_CLINICAL_SET_MT.out.vcf
                    .set { ch_clinical_mtsnv_vcf }
            }

            ch_ann_csq_mtsnv_in = ch_clinical_mtsnv_vcf.mix(ch_clin_research_mt_vcf.research)

            ANN_CSQ_PLI_MT(
                ch_variant_consequences_snv,
                ch_ann_csq_mtsnv_in,
                false
            )

            ANN_CSQ_PLI_MT.out.vcf_ann
                .filter { meta, _vcf ->
                    if (meta.probands.size()==0) {
                        log.warn("Skipping mitochondrial SNV ranking since no affected samples are detected in the case")
                    }
                    meta.probands.size()>0
                }
                .set {ch_ranksnv_mt_in}

            RANK_VARIANTS_MT (
                ch_pedfile,
                ch_reduced_penetrance,
                ch_score_config_mt,
                ch_ranksnv_mt_in,
                false
            )
        }
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL AND ANNOTATE NUCLEAR AND MITOCHONDRIAL SVs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    if (!skip_sv_calling) {
        CALL_STRUCTURAL_VARIANTS (
            ch_mapped.genome_marked_bam,
            ch_mapped.genome_marked_bai,
            ch_mapped.genome_marked_bam_bai,
            ch_genome_bwaindex,
            ch_genome_fasta,
            ch_genome_fai,
            ch_case_info,
            ch_target_bed,
            ch_genome_dictionary,
            ch_svcaller_priority,
            ch_readcount_intervals,
            ch_ploidy_model,
            ch_gcnvcaller_model,
            val_analysis_type,
            skip_germlinecnvcaller,
            ch_mapped.genome_marked_bam_bai,
            ch_genome_chrsizes,
            ch_genome_hisat2index,
            ch_mt_fai,
            ch_mt_fasta,
            ch_mt_lastdb,
            ch_input_fastqs,
            ch_subdepth,
            params.heavy_strand_origin_start,
            params.heavy_strand_origin_end,
            params.light_strand_origin_start,
            params.light_strand_origin_end,
            params.mito_length,
            params.mito_name,
            params.mitosalt_breakspan,
            params.mitosalt_breakthreshold,
            params.mitosalt_cluster_threshold,
            params.mitosalt_deletion_threshold_max,
            params.mitosalt_deletion_threshold_min,
            params.mitosalt_evalue_threshold,
            params.mitosalt_exclude,
            params.mitosalt_flank,
            params.mitosalt_heteroplasmy_limit,
            params.mitosalt_paired_distance,
            params.mitosalt_score_threshold,
            params.mitosalt_sizelimit,
            params.mitosalt_split_distance_threshold,
            params.mitosalt_split_length,
            params.saltshaker_dominant_fraction,
            params.saltshaker_group_radius,
            params.saltshaker_high_heteroplasmy,
            params.saltshaker_multiple_threshold,
            params.saltshaker_noise_threshold,
            val_run_mt_for_wes
        )
        ch_versions = ch_versions.mix(CALL_STRUCTURAL_VARIANTS.out.versions)

        //
        // ANNOTATE STRUCTURAL VARIANTS
        //
        if (!skip_sv_annotation) {
            ANNOTATE_STRUCTURAL_VARIANTS (
                ch_genome_dictionary,
                ch_genome_fasta,
                ch_svdb_bedpedbs,
                ch_svdb_dbs,
                CALL_STRUCTURAL_VARIANTS.out.vcf,
                ch_vep_cache,
                ch_vep_extra_files,
                val_svdb_query_bedpedbs,
                val_svdb_query_dbs,
                val_genome,
                val_vep_cache_version
            ).set { ch_sv_annotate }

            ch_sv_annotate.vcf_ann
                .multiMap { meta, vcf ->
                    clinical: [ meta + [ set: "clinical" ], vcf ]
                    research: [ meta + [ set: "research" ], vcf ]
                }
                .set { ch_clin_research_sv_vcf }

            ch_clinical_sv_vcf = channel.empty()
            if (!skip_generate_clinical_set) {
                GENERATE_CLINICAL_SET_SV(
                    ch_clin_research_sv_vcf.clinical,
                    ch_hgnc_ids,
                    false,
                    true
                )
                GENERATE_CLINICAL_SET_SV.out.vcf
                .set { ch_clinical_sv_vcf }
            }

            ch_ann_csq_sv_in = ch_clinical_sv_vcf.mix(ch_clin_research_sv_vcf.research)

            ANN_CSQ_PLI_SV (
                ch_variant_consequences_sv,
                ch_ann_csq_sv_in,
                false
            )

            ANN_CSQ_PLI_SV.out.vcf_ann
                .filter { meta, _vcf ->
                    if (meta.probands.size()==0) {
                        log.warn("Skipping SV ranking since no affected samples are detected in the case")
                    }
                    meta.probands.size()>0
                }
                .set {ch_ranksnv_sv_in}

            RANK_VARIANTS_SV (
                ch_pedfile,
                ch_reduced_penetrance,
                ch_score_config_sv,
                ch_ranksnv_sv_in,
                true
            )
        }
    }
/*

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL AND ANNOTATE MOBILE ELEMENTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    if (!skip_me_calling && val_analysis_type.equals("wgs")) {
        CALL_MOBILE_ELEMENTS(
            ch_case_info,
            ch_mapped.genome_marked_bam_bai,
            ch_genome_fai,
            ch_genome_fasta,
            ch_me_references
        )

        if (!skip_me_annotation) {
            ANNOTATE_MOBILE_ELEMENTS(
                ch_genome_dictionary,
                ch_genome_fasta,
                ch_me_svdb_resources,
                CALL_MOBILE_ELEMENTS.out.vcf,
                ch_vep_cache,
                val_genome,
                val_vep_cache_version,
                ch_vep_extra_files
            ).set { ch_me_annotate }

            ch_me_annotate.vcf_ann
                .multiMap { meta, vcf ->
                    clinical: [ meta + [ set: "clinical" ], vcf ]
                    research: [ meta + [ set: "research" ], vcf ]
                }
                .set { ch_clin_research_me_vcf }

            ch_clinical_me_vcf = channel.empty()
            if (!skip_generate_clinical_set) {
                GENERATE_CLINICAL_SET_ME(
                    ch_clin_research_me_vcf.clinical,
                    ch_hgnc_ids,
                    false,
                    true
                )
                GENERATE_CLINICAL_SET_ME.out.vcf
                .set { ch_clinical_me_vcf }
            }

            ch_ann_csq_me_in = ch_clinical_me_vcf.mix(ch_clin_research_me_vcf.research)

            ANN_CSQ_PLI_ME(
                ch_variant_consequences_sv,
                ch_ann_csq_me_in,
                true
            )

        }
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SMNCOPYNUMBERCALLER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    if ( val_analysis_type.equals("wgs") && !skip_smncopynumbercaller ) {

        RENAME_BAM.out.output
            .collect{_meta, bam -> bam}
            .toList()
            .set { ch_bam_list }

        RENAME_BAI.out.output
            .collect{_meta, bai -> bai}
            .toList()
            .set { ch_bai_list }

        ch_case_info
            .combine(ch_bam_list)
            .combine(ch_bai_list)
            .set { ch_bams_bais }

        SMNCOPYNUMBERCALLER (
            ch_bams_bais
        )
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PEDDY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    if (!skip_peddy) {
        PEDDY (
            CALL_SNV.out.genome_vcf.join(CALL_SNV.out.genome_tabix, failOnMismatch:true, failOnDuplicate:true),
            ch_pedfile.map{ped -> return[[id:"pedigree"], ped]},
            [[:],[]]
        )
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Generate CGH files from sequencing data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    if (!skip_vcf2cytosure && val_analysis_type.equals("wgs") && !skip_sv_calling && !skip_sv_annotation) {
        GENERATE_CYTOSURE_FILES (
            ch_mapped.genome_marked_bam_bai,
            ch_vcf2cytosure_blacklist,
            ch_sample_id_map,
            ch_sv_annotate.tbi,
            ch_sv_annotate.vcf_ann,
            val_sample_id_map
        )
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    if (!skip_gens && val_analysis_type.equals("wgs")) {
        GENS (
            ch_mapped.genome_marked_bam_bai,
            ch_genome_dictionary,
            ch_genome_fai,
            ch_genome_fasta,
            ch_gens_gnomad_pos,
            CALL_SNV.out.genome_gvcf,
            ch_gens_interval_list,
            ch_gens_pon_female,
            ch_gens_pon_male
        )
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VARIANT EVALUATION WITH RTGTOOLS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    if (val_run_rtgvcfeval) {
        VARIANT_EVALUATION (
            ch_rtg_truthvcfs,
            ch_sdf,
            CALL_SNV.out.genome_vcf_tabix
        )
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT SOFTWARE VERSIONS & MultiQC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'raredisease_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.fromPath("$projectDir/docs/images/nf-core-raredisease_logo_light.png", checkIfExists: true)


    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    if (!skip_fastqc) {
        ch_multiqc_files = ch_multiqc_files.mix(fastqc_report.collect{_meta, reports -> reports}.ifEmpty([]))
    }
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN.out.fastp_json.map{_meta, reports -> reports}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN.out.markdup_metrics.map{_meta, reports -> reports}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.sex_check.map{_meta, reports -> reports}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.multiple_metrics.map{_meta, reports -> reports}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.hs_metrics.map{_meta, reports -> reports}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.qualimap_results.map{_meta, reports -> reports}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.global_dist.map{_meta, reports -> reports}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.cov.map{_meta, reports -> reports}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.self_sm.map{_meta, reports -> reports}.collect().ifEmpty([]))

    if (!skip_peddy) {
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.ped.map{_meta, reports -> reports}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.het_check_csv.map{_meta, reports -> reports}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.ped_check_csv.map{_meta, reports -> reports}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.sex_check_csv.map{_meta, reports -> reports}.collect().ifEmpty([]))
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        ch_multiqc_samples
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
