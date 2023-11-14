//
// Convert VCF with structural variations to the “.CGH” format used by the CytoSure Interpret Software
//

include { TIDDIT_COV   } from '../../modules/nf-core/tiddit/cov/main'
include { VCF2CYTOSURE } from '../modules/nf-core/vcf2cytosure/main'

workflow SV_VCF_TO_CYTOSURE {
    take:
        ch_vcf  // channel: [mandatory] [ val(meta), path(vcf), path(vcf_index) ]
        ch_bam   // channel: [mandatory] [ path(bam) ]

    main:
        ch_versions = Channel.empty()

        TIDDIT_COV (ch_bam, [[],[]]) // 2nd pos. arg is req. only for cram input

        ch_vcf.dump(tag: 'vcf_channel', pretty: true)



        // Split vcf in samples combine and drop 

        // run vcf2cytosure samplevise


        [id:sample, [vcf]]


        ch_versions = ch_versions.mix(TIDDIT_COV.out.versions.first())

    emit:
        versions = ch_versions   // channel: [ versions.yml ]
}
