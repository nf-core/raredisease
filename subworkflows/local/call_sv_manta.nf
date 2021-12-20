//
// A structural variant caller workflow for manta
//

include { MANTA_GERMLINE as MANTA } from '../../modules/local/manta/germline/main'

workflow CALL_SV_MANTA {
    take:
    bam            // channel: [ val(meta), path(bam) ]
    bai            // channel: [ val(meta), path(bai) ]
    fasta          // path(fasta)
    fai            // path(fai)
    ch_case_info   // channel: [ case_id ]
    ch_bed         // channel: [ val(meta), path(bed), path(bed_tbi) ]

    main:
        bam.collect{it[1]}
            .toList()
            .set { bam_file_list }

        bai.collect{it[1]}
            .set { bai_file_list }

        ch_case_info.combine(bam_file_list)
            .set { manta_input_bams }

        ch_target_bed = ch_bed.ifEmpty([[],[],[]])

        MANTA ( manta_input_bams, bai_file_list, fasta, fai, ch_target_bed )
        ch_versions = MANTA.out.versions

    emit:
        candidate_small_indels_vcf      = MANTA.out.candidate_small_indels_vcf
        candidate_small_indels_vcf_tbi  = MANTA.out.candidate_small_indels_vcf_tbi
        candidate_sv_vcf                = MANTA.out.candidate_sv_vcf
        candidate_sv_vcf_tbi            = MANTA.out.candidate_sv_vcf_tbi
        diploid_sv_vcf                  = MANTA.out.diploid_sv_vcf
        diploid_sv_vcf_tbi              = MANTA.out.diploid_sv_vcf_tbi
        versions                        = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
