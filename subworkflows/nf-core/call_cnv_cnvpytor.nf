//
// CNVpytor workflow - Calling CNVs
//

include { CNVPYTOR_IMPORTREADDEPTH as GENERATE_PYTOR } from '../../modules/nf-core/cnvpytor/importreaddepth/main'
include { CNVPYTOR_HISTOGRAM as HISTOGRAMS           } from '../../modules/nf-core/cnvpytor/histogram/main'
include { CNVPYTOR_PARTITION as PARTITIONS           } from '../../modules/nf-core/cnvpytor/partition/main'
include { CNVPYTOR_CALLCNVS as CALL_CNVS             } from '../../modules/nf-core/cnvpytor/callcnvs/main'
include { CNVPYTOR_VIEW as VIEW                      } from '../../modules/nf-core/cnvpytor/view/main'

workflow CALL_CNV_CNVPYTOR {
    take:
        bam            // channel: [ val(meta), path(bam)]
        bai            // channel: [ val(meta), path(bai) ]
        case_info      // channel: [ case_id ]
        binsizes       // channel: [ val(binsize) ]
        fasta          // channel: [ path(fasta) ]
        fai            // channel: [ path(fai) ]


    main:
        ch_versions = Channel.empty()

        GENERATE_PYTOR(bam.join(bai, by: [0]), fasta, fai)
        ch_versions = ch_versions.mix(GENERATE_PYTOR.out.versions.first())

        HISTOGRAMS(GENERATE_PYTOR.out.pytor, binsizes)
        ch_versions = ch_versions.mix(HISTOGRAMS.out.versions.first())

        PARTITIONS(HISTOGRAMS.out.pytor, binsizes)
        ch_versions = ch_versions.mix(PARTITIONS.out.versions.first())

        CALL_CNVS(PARTITIONS.out.pytor, binsizes)
        ch_versions = ch_versions.mix(CALL_CNVS.out.versions.first())

        CALL_CNVS.out
            .pytor
            .collect{it[1]}
            .toList()
            .set { file_list }

        case_info
            .combine(file_list)
            .set { ch_pytor }

        VIEW(ch_pytor, binsizes, "vcf")
        ch_versions = ch_versions.mix(VIEW.out.versions.first())

    emit:
        candidate_cnvs_vcf              = VIEW.out.vcf              // channel: [ val(meta), path(*.tsv) ]
        versions                        = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}

