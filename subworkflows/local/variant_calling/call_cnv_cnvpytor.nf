//
// CNVpytor workflow - Calling CNVs
//

include { CNVPYTOR_IMPORTREADDEPTH as GENERATE_PYTOR } from '../../../modules/nf-core/cnvpytor/importreaddepth/main'
include { CNVPYTOR_HISTOGRAM as HISTOGRAMS           } from '../../../modules/nf-core/cnvpytor/histogram/main'
include { CNVPYTOR_PARTITION as PARTITIONS           } from '../../../modules/nf-core/cnvpytor/partition/main'
include { CNVPYTOR_CALLCNVS as CALL_CNVS             } from '../../../modules/nf-core/cnvpytor/callcnvs/main'
include { CNVPYTOR_VIEW as VIEW                      } from '../../../modules/nf-core/cnvpytor/view/main'

workflow CALL_CNV_CNVPYTOR {
    take:
        ch_bam            // channel: [mandatory] [ val(meta), path(bam)]
        ch_bai            // channel: [mandatory] [ val(meta), path(bai) ]
        ch_case_info      // channel: [mandatory] [ val(case_id) ]
        val_binsizes      // channel: [optional] [ val(binsize) ]
        ch_fasta          // channel: [mandatory] [ path(fasta) ]
        ch_fai            // channel: [mandatory] [ path(fai) ]


    main:
        ch_versions = Channel.empty()

        GENERATE_PYTOR(ch_bam.join(ch_bai, by: [0]), ch_fasta, ch_fai)

        HISTOGRAMS(GENERATE_PYTOR.out.pytor, binsizes)

        PARTITIONS(HISTOGRAMS.out.pytor, binsizes)

        CALL_CNVS(PARTITIONS.out.pytor, binsizes)

        CALL_CNVS.out
            .pytor
            .collect{it[1]}
            .toList()
            .set { file_list }

        ch_case_info
            .combine(file_list)
            .set { ch_pytor }

        VIEW(ch_pytor, val_binsizes, "vcf")

        ch_versions = ch_versions.mix(GENERATE_PYTOR.out.versions.first())
        ch_versions = ch_versions.mix(HISTOGRAMS.out.versions.first())
        ch_versions = ch_versions.mix(PARTITIONS.out.versions.first())
        ch_versions = ch_versions.mix(CALL_CNVS.out.versions.first())
        ch_versions = ch_versions.mix(VIEW.out.versions.first())

    emit:
        candidate_cnvs_vcf   = VIEW.out.vcf    // channel: [ val(meta), path(vcf) ]
        versions             = ch_versions     // channel: [ path(versions.yml) ]
}

