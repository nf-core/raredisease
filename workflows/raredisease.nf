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
include { PARSE_CONTAMINATION              } from '../modules/local/parse_contamination/main'
include { SANITY_CHECK_VCFANNO_DATABASES   } from '../modules/local/sanity_check_vcfanno_databases/main'

//
// SUBWORKFLOWS
//

include { ALIGN                                                       } from '../subworkflows/local/align'
include { ANNOTATE_CSQ_PLI as ANN_CSQ_PLI_ME                          } from '../subworkflows/local/annotate_consequence_pli'
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
include { CONTAMINATION_CHECK                                         } from '../subworkflows/local/contamination_check/main'
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
include { VCF_FILTER_BCFTOOLS_FILTERVEP as GENERATE_CLINICAL_SET_ME  } from '../subworkflows/local/vcf_filter_bcftools_filtervep'
include { VCF_FILTER_BCFTOOLS_FILTERVEP as GENERATE_CLINICAL_SET_MT  } from '../subworkflows/local/vcf_filter_bcftools_filtervep'
include { VCF_FILTER_BCFTOOLS_FILTERVEP as GENERATE_CLINICAL_SET_SNV } from '../subworkflows/local/vcf_filter_bcftools_filtervep'
include { VCF_FILTER_BCFTOOLS_FILTERVEP as GENERATE_CLINICAL_SET_SV  } from '../subworkflows/local/vcf_filter_bcftools_filtervep'

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
    ch_cadd_prescored
    ch_cadd_resources
    ch_call_interval
    ch_case_info
    ch_dbsnp
    ch_dbsnp_tbi
    ch_foundin_header
    ch_gcnvcaller_model
    ch_genome_bwafastalignindex
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
    ch_manta_regions
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
    ch_scatter_genome_split_intervals
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
    skip_mitosalt
    skip_ngsbits
    skip_peddy
    skip_smncopynumbercaller
    skip_vcf2cytosure
    val_aligner
    val_analysis_type
    val_cadd_resources
    val_concatenate_snv_calls
    val_skip_split_multiallelics
    val_exclude_alt
    val_extract_alignments
    val_genome
    val_heavy_strand_origin_end
    val_heavy_strand_origin_start
    val_homoplasmy_af_threshold
    val_light_strand_origin_end
    val_light_strand_origin_start
    val_mito_length
    val_mito_name
    val_mitosalt_breakspan
    val_mitosalt_breakthreshold
    val_mitosalt_cluster_threshold
    val_mitosalt_deletion_threshold_max
    val_mitosalt_deletion_threshold_min
    val_mitosalt_evalue_threshold
    val_mitosalt_exclude
    val_mitosalt_flank
    val_mitosalt_heteroplasmy_limit
    val_mitosalt_paired_distance
    val_mitosalt_score_threshold
    val_mitosalt_sizelimit
    val_mitosalt_split_distance_threshold
    val_mitosalt_split_length
    val_mt_aligner
    val_mt_subsample_approach
    val_mt_subsample_rd
    val_mt_subsample_seed
    val_multiqc_config
    val_multiqc_logo
    val_multiqc_methods_description
    val_multiqc_samples
    val_outdir
    val_platform
    val_run_mt_for_wes
    val_run_rtgvcfeval
    val_run_vcfanno_db_sanity_check
    val_sample_id_map
    val_save_all_mapped_as_cram
    val_save_noalt_mapped_as_cram
    val_svdb_query_bedpedbs
    val_svdb_query_dbs
    val_target_bed
    val_variant_caller
    val_vep_cache_version
    skip_contamination
    ch_contamination_sites
    ch_intervals_contamination

    main:

    ch_multiqc_files                    = channel.empty()
    ch_align_fastp_out                  = channel.empty()
    ch_align_genome_marked_cram         = channel.empty()
    ch_align_genome_marked_crai         = channel.empty()
    ch_align_genome_marked_bam          = channel.empty()
    ch_align_genome_marked_bai          = channel.empty()
    ch_align_markdup_metrics            = channel.empty()
    ch_subsample_mt_bam                 = channel.empty()
    ch_subsample_mt_bai                 = channel.empty()
    ch_call_snv_publish                 = channel.empty()
    ch_call_sv_vcf                      = channel.empty()
    ch_call_sv_tbi                      = channel.empty()
    ch_mt_del_result                    = channel.empty()
    ch_saltshaker_html                  = channel.empty()
    ch_saltshaker_plot                  = channel.empty()
    ch_call_snv_bcftools_concat_csi     = channel.empty()
    ch_call_snv_bcftools_concat_tbi     = channel.empty()
    ch_call_snv_bcftools_concat_vcf     = channel.empty()
    ch_call_snv_deepvariant_report      = channel.empty()
    ch_call_snv_genome_tabix            = channel.empty()
    ch_call_snv_genome_vcf              = channel.empty()
    ch_call_snv_mt_tabix                = channel.empty()
    ch_call_snv_mt_vcf                  = channel.empty()
    ch_call_sv_publish                  = channel.empty()
    ch_call_repeat_expansions_expansionhunter_bai = channel.empty()
    ch_call_repeat_expansions_expansionhunter_bam = channel.empty()
    ch_call_repeat_expansions_expansionhunter_vcf = channel.empty()
    ch_call_repeat_expansions_stranger_tbi        = channel.empty()
    ch_call_repeat_expansions_stranger_vcf        = channel.empty()
    ch_call_mobile_elements_tbi         = channel.empty()
    ch_call_mobile_elements_vcf         = channel.empty()
    ch_ann_csq_pli_me_tbi               = channel.empty()
    ch_ann_csq_pli_me_vcf_ann           = channel.empty()
    ch_annotate_genome_snvs_bcftools_concat_tbi       = channel.empty()
    ch_annotate_genome_snvs_bcftools_concat_vcf       = channel.empty()
    ch_annotate_genome_snvs_chromograph_autozyg_plots = channel.empty()
    ch_annotate_genome_snvs_chromograph_regions_plots = channel.empty()
    ch_annotate_genome_snvs_chromograph_sites_plots   = channel.empty()
    ch_annotate_genome_snvs_rhocall_viz_bed           = channel.empty()
    ch_annotate_genome_snvs_rhocall_viz_wig           = channel.empty()
    ch_annotate_genome_snvs_ucsc_wigtobigwig_bw       = channel.empty()
    ch_annotate_mt_snvs_ensemblvep_mt_tbi             = channel.empty()
    ch_annotate_mt_snvs_ensemblvep_mt_vcf             = channel.empty()
    ch_annotate_sv_publish              = channel.empty()
    ch_generate_cytosure_files_publish  = channel.empty()
    ch_gens_publish                     = channel.empty()
    ch_fastqc_publish                   = channel.empty()
    ch_smncopynumbercaller_publish      = channel.empty()
    ch_peddy_publish                    = channel.empty()
    ch_multiqc_publish                  = channel.empty()
    ch_rank_snv_publish                 = channel.empty()
    ch_rank_mt_publish                  = channel.empty()
    ch_rank_sv_publish                  = channel.empty()
    ch_variant_evaluation_publish       = channel.empty()

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
        ch_fastqc_publish = FASTQC.out.html
            .mix(FASTQC.out.zip)
            .map { meta, value -> ["fastqc/${meta.id}/", [meta, value]] }
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALIGN & FETCH STATS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ALIGN (
        ch_alignments,
        ch_genome_bwafastalignindex,
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
        val_exclude_alt,
        val_extract_alignments,
        val_mt_aligner,
        val_platform,
        val_run_mt_for_wes,
        val_save_all_mapped_as_cram,
        val_save_noalt_mapped_as_cram
    )
    .set { ch_mapped }
    ch_align_fastp_out              = ALIGN.out.fastp_out
    ch_align_genome_marked_bam      = ALIGN.out.genome_marked_bam
    ch_align_genome_marked_bai      = ALIGN.out.genome_marked_bai
    ch_align_genome_marked_cram     = ALIGN.out.genome_marked_cram
    ch_align_genome_marked_crai     = ALIGN.out.genome_marked_crai
    ch_align_markdup_metrics        = ALIGN.out.markdup_metrics

    if (!(skip_mt_subsample) && (val_analysis_type.equals("wgs") || val_run_mt_for_wes)) {
        if (val_mt_subsample_approach.equals("fraction")) {
            SUBSAMPLE_MT_FRAC(
                ch_mapped.mt_bam_bai,
                val_mt_subsample_rd,
                val_mt_subsample_seed
            )
            ch_subsample_mt_bam = SUBSAMPLE_MT_FRAC.out.bam
            ch_subsample_mt_bai = SUBSAMPLE_MT_FRAC.out.bai
        } else {
            SUBSAMPLE_MT_READS(
                ch_mapped.mt_bam_bai,
            )
            ch_subsample_mt_bam = SUBSAMPLE_MT_READS.out.bam
            ch_subsample_mt_bai = SUBSAMPLE_MT_READS.out.bai
        }
    }

    //
    // BAM QUALITY CHECK
    //
    QC_BAM (
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
        skip_ngsbits
    )

    //
    // SUBWORKFLOW: Check for contamination using GATK
    //
    ch_contamination_mqc    = Channel.empty()
    ch_contamination_table  = Channel.empty()
    ch_contamination_pileup = Channel.empty()

    if (!skip_contamination) {

        // Prepare BAM input with BAI
        ch_bam_for_contamination = ch_mapped.genome_marked_bam
            .join(ch_mapped.genome_marked_bai, failOnMismatch:true, failOnDuplicate:true)

        CONTAMINATION_CHECK (
            ch_bam_for_contamination,
            ch_genome_fasta,
            ch_genome_fai,
            ch_genome_dictionary,
            ch_contamination_sites,
            ch_intervals_contamination
        )

        // Parse for MultiQC
        PARSE_CONTAMINATION (
            CONTAMINATION_CHECK.out.contamination_table
        )

        ch_contamination_mqc    = PARSE_CONTAMINATION.out.mqc_table
        ch_contamination_table  = CONTAMINATION_CHECK.out.contamination_table
        ch_contamination_pileup = CONTAMINATION_CHECK.out.pileup_table
    }
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
        ch_call_repeat_expansions_expansionhunter_bai = CALL_REPEAT_EXPANSIONS.out.expansionhunter_bai
        ch_call_repeat_expansions_expansionhunter_bam = CALL_REPEAT_EXPANSIONS.out.expansionhunter_bam
        ch_call_repeat_expansions_expansionhunter_vcf = CALL_REPEAT_EXPANSIONS.out.expansionhunter_vcf

        if (!skip_repeat_annotation) {
            STRANGER (
                CALL_REPEAT_EXPANSIONS.out.vcf,
                ch_variant_catalog
            )
            ch_call_repeat_expansions_stranger_vcf = STRANGER.out.vcf
            ch_call_repeat_expansions_stranger_tbi = STRANGER.out.tbi
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
            val_skip_split_multiallelics,
            val_variant_caller
        )
        ch_call_snv_bcftools_concat_csi = CALL_SNV.out.bcftools_concat_csi
        ch_call_snv_bcftools_concat_tbi = CALL_SNV.out.bcftools_concat_tbi
        ch_call_snv_bcftools_concat_vcf = CALL_SNV.out.bcftools_concat_vcf
        ch_call_snv_deepvariant_report  = CALL_SNV.out.deepvariant_report
        ch_call_snv_genome_tabix        = CALL_SNV.out.genome_tabix
        ch_call_snv_genome_vcf          = CALL_SNV.out.genome_vcf
        ch_call_snv_mt_tabix            = CALL_SNV.out.mt_tabix
        ch_call_snv_mt_vcf              = CALL_SNV.out.mt_vcf

        // Removes vcfanno resource with empty records to keep vcfanno from crashing on those files
        ch_vcfanno_toml_final = ch_vcfanno_toml
        if (val_run_vcfanno_db_sanity_check && (!skip_snv_annotation || (!skip_mt_annotation && (val_run_mt_for_wes || val_analysis_type.matches("wgs|mito"))))) {
            ch_vcfanno_resources
                .combine(ch_vcfanno_extra)
                .map { files -> files.flatten() }
                .set { ch_all_vcfanno_dbs }
            SANITY_CHECK_VCFANNO_DATABASES (ch_vcfanno_toml, ch_all_vcfanno_dbs)
            ch_vcfanno_toml_final = SANITY_CHECK_VCFANNO_DATABASES.out.toml.collect()
        }

        //
        // ANNOTATE GENOME SNVs
        //
        if (!skip_snv_annotation) {

            ANNOTATE_GENOME_SNVS (
                ch_cadd_header,
                ch_cadd_prescored,
                ch_cadd_resources,
                ch_genome_chrsizes,
                ch_genome_fai,
                ch_genome_fasta,
                ch_gnomad_af,
                ch_samples,
                ch_scatter_genome_split_intervals,
                CALL_SNV.out.genome_vcf_tabix,
                ch_vcfanno_extra,
                ch_vcfanno_lua,
                ch_vcfanno_resources,
                ch_vcfanno_toml_final,
                ch_vep_cache,
                ch_vep_extra_files,
                val_analysis_type,
                val_cadd_resources,
                val_genome,
                val_vep_cache_version
            )
            ch_annotate_genome_snvs_bcftools_concat_tbi       = ANNOTATE_GENOME_SNVS.out.bcftools_concat_tbi
            ch_annotate_genome_snvs_bcftools_concat_vcf       = ANNOTATE_GENOME_SNVS.out.bcftools_concat_vcf
            ch_annotate_genome_snvs_chromograph_autozyg_plots = ANNOTATE_GENOME_SNVS.out.chromograph_autozyg_plots
            ch_annotate_genome_snvs_chromograph_regions_plots = ANNOTATE_GENOME_SNVS.out.chromograph_regions_plots
            ch_annotate_genome_snvs_chromograph_sites_plots   = ANNOTATE_GENOME_SNVS.out.chromograph_sites_plots
            ch_annotate_genome_snvs_rhocall_viz_bed           = ANNOTATE_GENOME_SNVS.out.rhocall_viz_bed
            ch_annotate_genome_snvs_rhocall_viz_wig           = ANNOTATE_GENOME_SNVS.out.rhocall_viz_wig
            ch_annotate_genome_snvs_ucsc_wigtobigwig_bw       = ANNOTATE_GENOME_SNVS.out.ucsc_wigtobigwig_bw

            ch_annotate_genome_snvs_bcftools_concat_vcf
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
            ch_rank_snv_publish = RANK_VARIANTS_SNV.out.publish
        }

        //
        // ANNOTATE MT SNVs
        //
        if (!(skip_mt_annotation) && (val_run_mt_for_wes || val_analysis_type.matches("wgs|mito"))) {

            ANNOTATE_MT_SNVS (
                ch_cadd_header,
                ch_cadd_prescored,
                ch_cadd_resources,
                ch_genome_fasta,
                ch_genome_fai,
                CALL_SNV.out.mt_vcf_tbi,
                ch_vcfanno_extra,
                ch_vcfanno_lua,
                ch_vcfanno_resources,
                ch_vcfanno_toml_final,
                ch_vep_cache,
                ch_vep_extra_files,
                val_cadd_resources,
                val_genome,
                val_homoplasmy_af_threshold,
                val_vep_cache_version
            ).set { ch_mt_annotate }
            ch_annotate_mt_snvs_ensemblvep_mt_tbi = ch_mt_annotate.ensemblvep_mt_tbi
            ch_annotate_mt_snvs_ensemblvep_mt_vcf = ch_mt_annotate.ensemblvep_mt_vcf

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
            ch_rank_mt_publish = RANK_VARIANTS_MT.out.publish
        }
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CALL AND ANNOTATE NUCLEAR AND MITOCHONDRIAL SVs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    if (!skip_sv_calling) {
        channel.of([val_mitosalt_breakspan,
            val_mitosalt_breakthreshold,
            val_mitosalt_cluster_threshold,
            val_mitosalt_deletion_threshold_max,
            val_mitosalt_deletion_threshold_min,
            val_mitosalt_evalue_threshold,
            val_mitosalt_exclude,
            val_mitosalt_paired_distance,
            val_mitosalt_score_threshold,
            val_mitosalt_sizelimit,
            val_mitosalt_split_distance_threshold,
            val_mitosalt_split_length])
            .set{ ch_mitosalt_config }

        CALL_STRUCTURAL_VARIANTS (
            ch_genome_bwaindex,
            ch_case_info,
            ch_gcnvcaller_model,
            ch_mapped.genome_marked_bai,
            ch_mapped.genome_marked_bam,
            ch_mapped.genome_marked_bam_bai,
            ch_genome_chrsizes,
            ch_genome_dictionary,
            ch_genome_fai,
            ch_genome_fasta,
            ch_genome_hisat2index,
            ch_manta_regions,
            ch_mitosalt_config,
            ch_mapped.mt_bam_bai,
            ch_mt_fai,
            ch_mt_fasta,
            ch_mt_lastdb,
            ch_ploidy_model,
            ch_readcount_intervals,
            ch_input_fastqs,
            ch_sample_id_map,
            ch_subdepth,
            ch_svcaller_priority,
            skip_germlinecnvcaller,
            skip_mitosalt,
            val_analysis_type,
            val_heavy_strand_origin_end,
            val_heavy_strand_origin_start,
            val_light_strand_origin_end,
            val_light_strand_origin_start,
            val_mito_length,
            val_mito_name,
            val_mitosalt_flank,
            val_mitosalt_heteroplasmy_limit,
            val_run_mt_for_wes
        )
        ch_call_sv_vcf = CALL_STRUCTURAL_VARIANTS.out.vcf
        ch_call_sv_tbi = CALL_STRUCTURAL_VARIANTS.out.tbi
        ch_saltshaker_html = CALL_STRUCTURAL_VARIANTS.out.saltshaker_html
        ch_saltshaker_plot = CALL_STRUCTURAL_VARIANTS.out.saltshaker_plot
        ch_mt_del_result = CALL_STRUCTURAL_VARIANTS.out.mt_del_result

        //
        // ANNOTATE STRUCTURAL VARIANTS
        //
        if (!skip_sv_annotation) {
            ANNOTATE_STRUCTURAL_VARIANTS (
                ch_genome_dictionary,
                ch_genome_fasta,
                ch_svdb_bedpedbs,
                ch_svdb_dbs,
                ch_call_sv_vcf,
                ch_vep_cache,
                ch_vep_extra_files,
                val_svdb_query_bedpedbs,
                val_svdb_query_dbs,
                val_genome,
                val_vep_cache_version
            ).set { ch_sv_annotate }
            ch_annotate_sv_publish = ch_sv_annotate.publish

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
            ch_rank_sv_publish = RANK_VARIANTS_SV.out.publish
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
        ch_call_mobile_elements_vcf = CALL_MOBILE_ELEMENTS.out.vcf
        ch_call_mobile_elements_tbi = CALL_MOBILE_ELEMENTS.out.tbi

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
            ch_ann_csq_pli_me_vcf_ann = ANN_CSQ_PLI_ME.out.vcf_ann
            ch_ann_csq_pli_me_tbi     = ANN_CSQ_PLI_ME.out.tbi

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
        ch_smncopynumbercaller_publish = SMNCOPYNUMBERCALLER.out.smncopynumber
            .mix(SMNCOPYNUMBERCALLER.out.run_metrics)
            .map { meta, value -> ['smncopynumbercaller/', [meta, value]] }
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
        ch_peddy_publish = PEDDY.out.vs_html
            .mix(PEDDY.out.html)
            .mix(PEDDY.out.ped)
            .mix(PEDDY.out.het_check_png)
            .mix(PEDDY.out.ped_check_png)
            .mix(PEDDY.out.sex_check_png)
            .mix(PEDDY.out.het_check_csv)
            .mix(PEDDY.out.ped_check_csv)
            .mix(PEDDY.out.sex_check_csv)
            .mix(PEDDY.out.ped_check_rel_difference_csv)
            .map { meta, value -> ['peddy/', [meta, value]] }
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
        ch_generate_cytosure_files_publish = GENERATE_CYTOSURE_FILES.out.publish
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
            CALL_SNV.out.genome_gtabix,
            ch_gens_interval_list,
            ch_gens_pon_female,
            ch_gens_pon_male
        )
        ch_gens_publish = GENS.out.publish
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
        ch_variant_evaluation_publish = VARIANT_EVALUATION.out.publish
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

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${val_outdir}/pipeline_info",
            name: 'nf_core_'  +  'raredisease_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = val_multiqc_config ?
        channel.fromPath(val_multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = val_multiqc_logo ?
        channel.fromPath(val_multiqc_logo, checkIfExists: true) :
        channel.fromPath("$projectDir/docs/images/nf-core-raredisease_logo_light.png", checkIfExists: true)

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = val_multiqc_methods_description ?
        file(val_multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))
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
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.ngsbits_samplegender_tsv.map{_meta, reports -> reports}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.picard_collectmultiplemetrics_metrics.map{_meta, reports -> reports}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.picard_collecthsmetrics_metrics.map{_meta, reports -> reports}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.mosdepth_global_txt.map{_meta, reports -> reports}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.wgsmetrics_wg.map{_meta, reports -> reports}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_BAM.out.verifybamid_self_sm.map{_meta, reports -> reports}.collect().ifEmpty([]))

    // Add contamination results to MultiQC
    ch_multiqc_files = ch_multiqc_files.mix(ch_contamination_mqc.map { _meta, file -> file })

    if (!skip_peddy) {
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.ped.map{_meta, reports -> reports}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.het_check_csv.map{_meta, reports -> reports}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.ped_check_csv.map{_meta, reports -> reports}.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PEDDY.out.sex_check_csv.map{_meta, reports -> reports}.collect().ifEmpty([]))
    }

    // config: custom config if provided, otherwise the default pipeline config
    // logo: custom logo if provided, otherwise omitted
    // sample_names: optional TSV for MultiQC --sample-names
    MULTIQC (
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'raredisease'],
                files,
                val_multiqc_config
                    ? file(val_multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                val_multiqc_logo ? file(val_multiqc_logo, checkIfExists: true) : [],
                [],
                val_multiqc_samples ? file(val_multiqc_samples) : [],
            ]
        }
    )
    ch_multiqc_publish = MULTIQC.out.report
        .mix(MULTIQC.out.data)
        .mix(MULTIQC.out.plots)
        .map { _meta, value -> ['multiqc/', value] }

    emit:
    align_fastp_out              = ch_align_fastp_out              // channel: [ val(meta), path(json|html|log|reads|reads_fail|reads_merged) ]
    align_genome_marked_bam      = ch_align_genome_marked_bam      // channel: [ val(meta), path(bam) ]
    align_genome_marked_bai      = ch_align_genome_marked_bai      // channel: [ val(meta), path(bai) ]
    align_genome_marked_cram     = ch_align_genome_marked_cram     // channel: [ val(meta), path(cram) ]
    align_genome_marked_crai     = ch_align_genome_marked_crai     // channel: [ val(meta), path(crai) ]
    align_markdup_metrics        = ch_align_markdup_metrics        // channel: [ val(meta), path(metrics) ]
    multiqc_report                                   = MULTIQC.out.report.map { _meta, report -> report }.toList()
    qc_bam_chromograph_cov_plots                     = QC_BAM.out.chromograph_cov_plots                    // channel: [ val(meta), path(png) ]
    qc_bam_mosdepth_global_txt                       = QC_BAM.out.mosdepth_global_txt                      // channel: [ val(meta), path(txt) ]
    qc_bam_mosdepth_per_base_bed                     = QC_BAM.out.mosdepth_per_base_bed                    // channel: [ val(meta), path(bed.gz) ]
    qc_bam_mosdepth_per_base_csi                     = QC_BAM.out.mosdepth_per_base_csi                    // channel: [ val(meta), path(csi) ]
    qc_bam_mosdepth_per_base_d4                      = QC_BAM.out.mosdepth_per_base_d4                     // channel: [ val(meta), path(d4) ]
    qc_bam_mosdepth_quantized_bed                    = QC_BAM.out.mosdepth_quantized_bed                   // channel: [ val(meta), path(bed.gz) ]
    qc_bam_mosdepth_quantized_csi                    = QC_BAM.out.mosdepth_quantized_csi                   // channel: [ val(meta), path(csi) ]
    qc_bam_mosdepth_regions_bed                      = QC_BAM.out.mosdepth_regions_bed                     // channel: [ val(meta), path(bed.gz) ]
    qc_bam_mosdepth_regions_csi                      = QC_BAM.out.mosdepth_regions_csi                     // channel: [ val(meta), path(csi) ]
    qc_bam_mosdepth_regions_txt                      = QC_BAM.out.mosdepth_regions_txt                     // channel: [ val(meta), path(txt) ]
    qc_bam_mosdepth_summary_txt                      = QC_BAM.out.mosdepth_summary_txt                     // channel: [ val(meta), path(txt) ]
    qc_bam_mosdepth_thresholds_bed                   = QC_BAM.out.mosdepth_thresholds_bed                  // channel: [ val(meta), path(bed.gz) ]
    qc_bam_mosdepth_thresholds_csi                   = QC_BAM.out.mosdepth_thresholds_csi                  // channel: [ val(meta), path(csi) ]
    qc_bam_ngsbits_samplegender_tsv                  = QC_BAM.out.ngsbits_samplegender_tsv                 // channel: [ val(meta), path(tsv) ]
    qc_bam_picard_collecthsmetrics_metrics           = QC_BAM.out.picard_collecthsmetrics_metrics          // channel: [ val(meta), path(metrics) ]
    qc_bam_picard_collectmultiplemetrics_metrics     = QC_BAM.out.picard_collectmultiplemetrics_metrics    // channel: [ val(meta), path(metrics) ]
    qc_bam_picard_collectmultiplemetrics_pdf         = QC_BAM.out.picard_collectmultiplemetrics_pdf        // channel: [ val(meta), path(pdf) ]
    qc_bam_sambamba_depth_bed                        = QC_BAM.out.sambamba_depth_bed                       // channel: [ val(meta), path(bed) ]
    qc_bam_tiddit_cov_cov                            = QC_BAM.out.tiddit_cov_cov                           // channel: [ val(meta), path(bed) ]
    qc_bam_tiddit_cov_wig                            = QC_BAM.out.tiddit_cov_wig                           // channel: [ val(meta), path(wig) ]
    qc_bam_ucsc_wigtobigwig_bw                       = QC_BAM.out.ucsc_wigtobigwig_bw                      // channel: [ val(meta), path(bw) ]
    qc_bam_verifybamid_ancestry                      = QC_BAM.out.verifybamid_ancestry                     // channel: [ val(meta), path(ancestry) ]
    qc_bam_verifybamid_bed                           = QC_BAM.out.verifybamid_bed                          // channel: [ val(meta), path(bed) ]
    qc_bam_verifybamid_log                           = QC_BAM.out.verifybamid_log                          // channel: [ val(meta), path(log) ]
    qc_bam_verifybamid_mu                            = QC_BAM.out.verifybamid_mu                           // channel: [ val(meta), path(mu) ]
    qc_bam_verifybamid_self_sm                       = QC_BAM.out.verifybamid_self_sm                      // channel: [ val(meta), path(selfSM) ]
    qc_bam_verifybamid_ud                            = QC_BAM.out.verifybamid_ud                           // channel: [ val(meta), path(ud) ]
    qc_bam_wgsmetrics_wg                             = QC_BAM.out.wgsmetrics_wg                            // channel: [ val(meta), path(metrics) ]
    qc_bam_wgsmetrics_y                              = QC_BAM.out.wgsmetrics_y                             // channel: [ val(meta), path(metrics) ]
    call_sv_vcf                                      = ch_call_sv_vcf                                      // channel: [ val(meta), path(vcf) ]
    call_sv_tbi                                      = ch_call_sv_tbi                                      // channel: [ val(meta), path(tbi) ]
    saltshaker_html                                  = ch_saltshaker_html                                  // channel: [ val(meta), path(html) ]
    saltshaker_plot                                  = ch_saltshaker_plot                                  // channel: [ val(meta), path(png) ]
    mt_del_result                                    = ch_mt_del_result                                    // channel: [ val(meta), path(txt) ]
    call_repeat_expansions_expansionhunter_bai       = ch_call_repeat_expansions_expansionhunter_bai       // channel: [ val(meta), path(bai) ]
    call_repeat_expansions_expansionhunter_bam       = ch_call_repeat_expansions_expansionhunter_bam       // channel: [ val(meta), path(bam) ]
    call_repeat_expansions_expansionhunter_vcf       = ch_call_repeat_expansions_expansionhunter_vcf       // channel: [ val(meta), path(vcf) ]
    call_repeat_expansions_stranger_tbi              = ch_call_repeat_expansions_stranger_tbi              // channel: [ val(meta), path(tbi) ]
    call_repeat_expansions_stranger_vcf              = ch_call_repeat_expansions_stranger_vcf              // channel: [ val(meta), path(vcf) ]
    call_snv_bcftools_concat_csi             = ch_call_snv_bcftools_concat_csi                             // channel: [ val(meta), path(csi) ]
    call_snv_bcftools_concat_tbi             = ch_call_snv_bcftools_concat_tbi                             // channel: [ val(meta), path(tbi) ]
    call_snv_bcftools_concat_vcf             = ch_call_snv_bcftools_concat_vcf                             // channel: [ val(meta), path(vcf) ]
    call_snv_deepvariant_report              = ch_call_snv_deepvariant_report                              // channel: [ val(meta), path(html) ]
    call_snv_genome_tabix                    = ch_call_snv_genome_tabix                                    // channel: [ val(meta), path(tbi) ]
    call_snv_genome_vcf                      = ch_call_snv_genome_vcf                                      // channel: [ val(meta), path(vcf) ]
    call_snv_mt_tabix                        = ch_call_snv_mt_tabix                                        // channel: [ val(meta), path(tbi) ]
    call_snv_mt_vcf                          = ch_call_snv_mt_vcf                                          // channel: [ val(meta), path(vcf) ]
    annotate_genome_snvs_bcftools_concat_tbi         = ch_annotate_genome_snvs_bcftools_concat_tbi         // channel: [ val(meta), path(tbi) ]
    annotate_genome_snvs_bcftools_concat_vcf         = ch_annotate_genome_snvs_bcftools_concat_vcf         // channel: [ val(meta), path(vcf) ]
    annotate_genome_snvs_chromograph_autozyg_plots   = ch_annotate_genome_snvs_chromograph_autozyg_plots   // channel: [ val(meta), path(png) ]
    annotate_genome_snvs_chromograph_regions_plots   = ch_annotate_genome_snvs_chromograph_regions_plots   // channel: [ val(meta), path(png) ]
    annotate_genome_snvs_chromograph_sites_plots     = ch_annotate_genome_snvs_chromograph_sites_plots     // channel: [ val(meta), path(png) ]
    annotate_genome_snvs_rhocall_viz_bed             = ch_annotate_genome_snvs_rhocall_viz_bed             // channel: [ val(meta), path(bed) ]
    annotate_genome_snvs_rhocall_viz_wig             = ch_annotate_genome_snvs_rhocall_viz_wig             // channel: [ val(meta), path(wig) ]
    annotate_genome_snvs_ucsc_wigtobigwig_bw         = ch_annotate_genome_snvs_ucsc_wigtobigwig_bw         // channel: [ val(meta), path(bw) ]
    annotate_mt_snvs_ensemblvep_mt_tbi               = ch_annotate_mt_snvs_ensemblvep_mt_tbi // channel: [ val(meta), path(tbi) ]
    annotate_mt_snvs_ensemblvep_mt_vcf               = ch_annotate_mt_snvs_ensemblvep_mt_vcf // channel: [ val(meta), path(vcf) ]
    call_mobile_elements_tbi                         = ch_call_mobile_elements_tbi // channel: [ val(meta), path(tbi) ]
    call_mobile_elements_vcf                         = ch_call_mobile_elements_vcf // channel: [ val(meta), path(vcf) ]
    ann_csq_pli_me_tbi                               = ch_ann_csq_pli_me_tbi       // channel: [ val(meta), path(tbi) ]
    ann_csq_pli_me_vcf_ann                           = ch_ann_csq_pli_me_vcf_ann   // channel: [ val(meta), path(vcf) ]
    subsample_mt_bai             = ch_subsample_mt_bai             // channel: [ val(meta), path(bai) ]
    subsample_mt_bam             = ch_subsample_mt_bam             // channel: [ val(meta), path(bam) ]
    contamination_table          = ch_contamination_table         // channel: [ val(meta), path(table) ]
    contamination_pileup         = ch_contamination_pileup        // channel: [ val(meta), path(table) ]
    versions                     = ch_versions
    publish                      = ch_call_sv_publish
                       .mix(ch_annotate_sv_publish)
                       .mix(ch_generate_cytosure_files_publish)
                       .mix(ch_gens_publish)
                       .mix(ch_fastqc_publish)
                       .mix(ch_smncopynumbercaller_publish)
                       .mix(ch_peddy_publish)
                       .mix(ch_multiqc_publish)
                       .mix(ch_rank_snv_publish)
                       .mix(ch_rank_mt_publish)
                       .mix(ch_rank_sv_publish)
                       .mix(ch_variant_evaluation_publish)
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
