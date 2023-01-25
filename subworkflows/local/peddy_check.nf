//
// Peddy subworkflow to check sex and relatedness.
//

include { PEDDY } from '../modules/nf-core/peddy/main'

workflow PEDDY_CHECK {
    take:
        vcf                // channel: [ val(meta), path(vcf), path(vcf_index) ]
        ped                
        
    main:
        ch_versions = Channel.empty()

        PEDDY(vcf, ped)
        ch_versions = ch_versions.mix(PEDDY.out.versions)
        
    emit:
        ped                    = PEDDY.out.ped
        csv                    = PEDDY.out.csv
        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}


