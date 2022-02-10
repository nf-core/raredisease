//
// A nested subworkflow to call structural variants.
//

include { CALL_SV_MANTA } from './call_sv_manta'

include { CALL_SV_TIDDIT } from './call_sv_tiddit'

include { SVDB_MERGE } from '../../modules/nf-core/modules/svdb/merge/main'

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

        //manta
        CALL_SV_MANTA ( bam, bai, fasta, fai, case_info, target_bed ).diploid_sv_vcf
            .collect{it[1]}
            .set{ manta_vcf }
        ch_versions = ch_versions.mix(CALL_SV_MANTA.out.versions)

        //tiddit
        CALL_SV_TIDDIT ( bam, fasta, fai, case_info ).vcf
            .collect{it[1]}
            .set { tiddit_vcf }
        ch_versions = ch_versions.mix(CALL_SV_TIDDIT.out.versions)

        //merge
        tiddit_vcf
            .combine(manta_vcf)
            .toList()
            .set { vcf_list }

        case_info.combine(vcf_list)
            .set { merge_input_vcfs }

        SVDB_MERGE ( merge_input_vcfs, ["tiddit","manta"] )


    emit:
        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
