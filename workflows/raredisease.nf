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

    //
    // Prepare references and indices.
    //
    PREPARE_REFERENCES (
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_dictionary
    )
    .set { ch_references }

    //
    // Gather built indices or get them from the params
    //
    ch_genome_bwaindex          = params.bwa                                ? Channel.fromPath(params.bwa).map {it -> [[id:it.simpleName], it]}.collect()
                                                                            : ch_references.genome_bwa_index
    ch_genome_bwamem2index      = params.bwamem2                            ? Channel.fromPath(params.bwamem2).map {it -> [[id:it.simpleName], it]}.collect()
                                                                            : ch_references.genome_bwamem2_index
    ch_genome_bwamemeindex      = params.bwameme                            ? Channel.fromPath(params.bwameme).map {it -> [[id:it.simpleName], it]}.collect()
                                                                            : ch_references.genome_bwameme_index
    ch_genome_chrsizes          = ch_references.genome_chrom_sizes
    ch_genome_fai               = ch_references.genome_fai
    ch_genome_dictionary        = ch_references.genome_dict
    ch_versions                 = ch_versions.mix(ch_references.versions)

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
        params.mbuffer_mem,
        params.platform,
        params.samtools_sort_threads
    )
    .set { ch_mapped }
    ch_versions   = ch_versions.mix(ALIGN.out.versions)



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

    ch_multiqc_files = ch_multiqc_files.mix(ALIGN.out.fastp_json.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN.out.markdup_metrics.map{it[1]}.collect().ifEmpty([]))

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
