#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_PREPROCESSINTERVALS           } from '../../../modules/nf-core/gatk4/preprocessintervals/main.nf'
include { GATK4_ANNOTATEINTERVALS             } from '../../../modules/nf-core/gatk4/annotateintervals/main.nf'
include { GATK4_FILTERINTERVALS               } from '../../../modules/nf-core/gatk4/filterintervals/main.nf'
include { GATK4_INTERVALLISTTOOLS             } from '../../../modules/nf-core/gatk4/intervallisttools/main.nf'
include { GATK4_COLLECTREADCOUNTS             } from '../../../modules/nf-core/gatk4/collectreadcounts/main.nf'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY } from '../../../modules/nf-core/gatk4/determinegermlinecontigploidy/main.nf'
include { GATK4_GERMLINECNVCALLER             } from '../../../modules/nf-core/gatk4/germlinecnvcaller/main.nf'
include { GATK4_POSTPROCESSGERMLINECNVCALLS   } from '../../../modules/nf-core/gatk4/postprocessgermlinecnvcalls/main.nf'

workflow CALL_CNV_GERMLINECNVCALLER {
    take:
        bam_bai        // channel: [ val(meta), path(bam), path(bai) ]
        bam_bai_2      // channel: [ val(meta), path(bam), path(bai) ]
        fasta_no_meta  // channel: [ path(fasta_no_meta) ]
        fai            // channel: [ path(fai) ]
        target_bed     // channel: [ path(target_bed) ]
        blacklist_bed  // channel: [ val(meta), path(blacklist_bed) ]
        dict           // channel: [ path(dict) ]
        priors         // [ path(priors) ]
        ploidy_model   // channel: [ path(ploidy_model) ]
        cnv_model      // channel: [ path(cnv_model) ]

    main:
        ch_versions = Channel.empty()

        GATK4_PREPROCESSINTERVALS ( blacklist_bed, fasta_no_meta, fai, dict )

        GATK4_ANNOTATEINTERVALS ( GATK4_PREPROCESSINTERVALS.out.interval_list, fasta_no_meta, fai, dict, [], [], [], [])

        inputs = bam_bai.combine( GATK4_PREPROCESSINTERVALS.out.interval_list ).mix( bam_bai_2.combine( GATK4_PREPROCESSINTERVALS.out.interval_list ) )

        GATK4_COLLECTREADCOUNTS ( inputs, fasta_no_meta, fai, dict )

        GATK4_FILTERINTERVALS ( GATK4_PREPROCESSINTERVALS.out.interval_list, GATK4_COLLECTREADCOUNTS.out.tsv.collect{ it[1] }, GATK4_ANNOTATEINTERVALS.out.annotated_intervals )

        GATK4_INTERVALLISTTOOLS ( GATK4_FILTERINTERVALS.out.interval_list )

        dgcp_case_input = GATK4_COLLECTREADCOUNTS.out.tsv
                .map({ meta, tsv -> [ [id:'test'], tsv, [], [] ] })
                .groupTuple()
        GATK4_DETERMINEGERMLINECONTIGPLOIDY ( dgcp_case_input, [], ploidy_model )

        gcnvc_case_input = GATK4_COLLECTREADCOUNTS.out.tsv
                .map({ meta, tsv -> return [[id:"test"], tsv, [] ]})
                .groupTuple()
        GATK4_GERMLINECNVCALLER ( gcnvc_case_input, cnv_model, GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls.collect{ it[1] } )

        GATK4_POSTPROCESSGERMLINECNVCALLS ( GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls, cnv_model, GATK4_GERMLINECNVCALLER.out.calls.collect{ it[1] } )

        ch_versions = ch_versions.mix(GATK4_PREPROCESSINTERVALS.out.versions)
        ch_versions = ch_versions.mix(GATK4_ANNOTATEINTERVALS.out.versions)
        ch_versions = ch_versions.mix(GATK4_COLLECTREADCOUNTS.out.versions)
        ch_versions = ch_versions.mix(GATK4_FILTERINTERVALS.out.versions)
        ch_versions = ch_versions.mix(GATK4_INTERVALLISTTOOLS.out.versions)
        ch_versions = ch_versions.mix(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.versions)
        ch_versions = ch_versions.mix(GATK4_GERMLINECNVCALLER.out.versions)
        ch_versions = ch_versions.mix(GATK4_POSTPROCESSGERMLINECNVCALLS.out.versions)

    emit:
        candidate_cnvs_vcf_tar_gz       = GATK4_GERMLINECNVCALLER.out.tar_gz  // channel: [ val(meta), path(*.tar_gz) ]
        versions                        = ch_versions                         // channel: [ versions.yml ]
}
