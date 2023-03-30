//
// Peddy subworkflow to check sex and relatedness.
//

include { PEDDY } from '../../modules/nf-core/peddy/main'

workflow PEDDY_CHECK {
    take:
        ch_vcf  // channel: [mandatory] [ val(meta), path(vcf), path(vcf_index) ]
        ch_ped  // channel: [mandatory] [ path(ped) ]

    main:
        ch_versions = Channel.empty()

        PEDDY( ch_vcf, ch_ped )
        ch_versions = ch_versions.mix(PEDDY.out.versions.first())

    emit:
        ped      = PEDDY.out.ped // channel: [ val(meta), path(ped) ]
        csv      = PEDDY.out.csv // channel: [ val(meta), path(csv) ]
        versions = ch_versions   // channel: [ versions.yml ]
}
