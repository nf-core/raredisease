//
// A subworkflow to call CNVs using cnvnator
//

nextflow.enable.dsl = 2

include { CNVNATOR_CNVNATOR as CNVNATOR_RD        } from '../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_HIST      } from '../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_STAT      } from '../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_PARTITION } from '../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_CALL      } from '../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CONVERT2VCF                    } from '../../../modules/nf-core/cnvnator/convert2vcf/main.nf'
include { SVDB_MERGE as SVDB_MERGE_CNVNATOR       } from '../../../modules/nf-core/svdb/merge/main'

workflow CALL_SV_CNVNATOR {
    take:
        ch_bam_bai   // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_fasta     // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_fai       // channel: [mandatory] [ val(meta), path(fai) ]
        ch_case_info // channel: [mandatory] [ val(case_info) ]

    main:
        ch_versions = Channel.empty()

        CNVNATOR_RD ( ch_bam_bai, [[:],[]], [[:],[]], [[:],[]], "rd" )
        CNVNATOR_HIST ( [[:],[],[]], CNVNATOR_RD.out.root, ch_fasta, ch_fai, "his" )
        CNVNATOR_STAT ( [[:],[],[]], CNVNATOR_HIST.out.root, [[:],[]], [[:],[]], "stat" )
        CNVNATOR_PARTITION ( [[:],[],[]], CNVNATOR_STAT.out.root, [[:],[]], [[:],[]], "partition" )
        CNVNATOR_CALL ( [[:],[],[]], CNVNATOR_PARTITION.out.root, [[:],[]], [[:],[]], "call" )
        CNVNATOR_CONVERT2VCF (CNVNATOR_CALL.out.tab).vcf
            .collect{it[1]}
            .toList()
            .set { vcf_file_list }

        ch_case_info
            .combine(vcf_file_list)
            .set { merge_input_vcfs }

        SVDB_MERGE_CNVNATOR ( merge_input_vcfs, [] )

        ch_versions = ch_versions.mix(CNVNATOR_RD.out.versions)
        ch_versions = ch_versions.mix(CNVNATOR_HIST.out.versions)
        ch_versions = ch_versions.mix(CNVNATOR_STAT.out.versions)
        ch_versions = ch_versions.mix(CNVNATOR_PARTITION.out.versions)
        ch_versions = ch_versions.mix(CNVNATOR_CALL.out.versions)
        ch_versions = ch_versions.mix(CNVNATOR_CONVERT2VCF.out.versions)

    emit:
        vcf        = SVDB_MERGE_CNVNATOR.out.vcf  // channel: [ val(meta), path(*.tar.gz) ]
        versions   = ch_versions                  // channel: [ versions.yml ]
}
