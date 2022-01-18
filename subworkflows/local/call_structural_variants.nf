//
// A nested subworkflow to call structural variants.
//

include { CALL_SV_MANTA } from './call_sv_manta'

workflow CALL_STRUCTURAL_VARIANTS {

    take:
        bam         // channel: [ val(meta), path(bam) ]
        bai         // channel: [ val(meta), path(bai) ]
        fasta       // channel: [ path(genome.fasta) ]
        fai         // channel: [ path(genome.fai) ]
        case_info   // channel: [ val(case_info) ]
        target_bed  // channel: [ path(target.bed) ]

    main:
        ch_versions = Channel.empty()

        CALL_SV_MANTA( bam, bai, fasta, fai, case_info, target_bed )
        ch_versions = ch_versions.mix(CALL_SV_MANTA.out.versions)

    emit:
        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
