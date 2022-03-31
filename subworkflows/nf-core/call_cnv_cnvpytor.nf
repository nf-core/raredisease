//
// CNVpytor workflow - Calling CNVs
//

include { CNVPYTOR_IMPORTREADDEPTH as GENERATE_PYTOR } from '../../modules/nf-core/modules/cnvpytor/importreaddepth/main'
include { CNVPYTOR_HISTOGRAM as HISTOGRAMS           } from '../../modules/nf-core/modules/cnvpytor/histogram/main'
include { CNVPYTOR_PARTITION as PARTITIONS           } from '../../modules/nf-core/modules/cnvpytor/partition/main'
include { CNVPYTOR_CALLCNVS as CALL_CNVS             } from '../../modules/nf-core/modules/cnvpytor/callcnvs/main'

workflow CALL_CNV_CNVPYTOR {
    take:
        bam            // channel: [ val(meta), path(bam)]
        bai            // channel: [ val(meta), path(bai) ]
        case_info      // channel: [ case_id ]


    main:
        ch_versions = Channel.empty()

        GENERATE_PYTOR(bam.join(bai, by: [0]),[],[])
        HISTOGRAMS(GENERATE_PYTOR.out.pytor)
        PARTITIONS(HISTOGRAMS.out.pytor)
        CALL_CNVS(PARTITIONS.out.pytor)
        ch_versions = ch_versions.mix(CALL_CNVS.out.versions)
        //TO DO : tsv2vcf

    emit:
        candidate_cnvs_tsv              = CALL_CNVS.out.cnvs        // channel: [ val(meta), path(*.tsv) ]
        versions                        = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}

