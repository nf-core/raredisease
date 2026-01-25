//
// A structural variant caller workflow for manta
//

include { MANTA_GERMLINE as MANTA              } from '../../../modules/nf-core/manta/germline/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_MANTA } from '../../../modules/nf-core/bcftools/view/main.nf'

workflow CALL_SV_MANTA {
    take:
        ch_bam            // channel: [mandatory] [ val(meta), path(bam) ]
        ch_bai            // channel: [mandatory] [ val(meta), path(bai) ]
        ch_genome_fasta   // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai     // channel: [mandatory] [ val(meta), path(fai) ]
        ch_case_info      // channel: [mandatory] [ val(case_info) ]
        ch_bed            // channel: [mandatory for WES] [ val(meta), path(bed), path(tbi) ]
        val_analysis_type // string: "wes", "wgs", or "mito"

    main:
        ch_bam.collect{_meta, bam -> bam}
            .toList()
            .set { bam_file_list }

        ch_bai.collect{_meta, bai -> bai}
            .toList()
            .set { bai_file_list }

        ch_bed.map {
                _id, bed_file, index ->
                    return [bed_file, index]}
            .set { bed_input }

        if (val_analysis_type.equals("wgs")) {
            ch_case_info.combine(bam_file_list)
                .combine(bai_file_list)
                .map { meta, input, index -> [meta, input, index] + [ [], [] ] }
                .set { manta_input }
            MANTA ( manta_input, ch_genome_fasta, ch_genome_fai, [] )
        } else {
            ch_case_info.combine(bam_file_list)
                .combine(bai_file_list)
                .combine(bed_input)
                .set { manta_input }
            MANTA ( manta_input, ch_genome_fasta, ch_genome_fai, [] )
        }

        MANTA.out.diploid_sv_vcf
            .join(MANTA.out.diploid_sv_vcf_tbi)
            .set {ch_filter_in}
        BCFTOOLS_VIEW_MANTA (ch_filter_in, [], [], [])

        ch_versions = MANTA.out.versions

    emit:
        candidate_small_indels_vcf     = MANTA.out.candidate_small_indels_vcf     // channel: [ val(meta), path(vcf) ]
        candidate_small_indels_vcf_tbi = MANTA.out.candidate_small_indels_vcf_tbi // channel: [ val(meta), path(tbi) ]
        candidate_sv_vcf               = MANTA.out.candidate_sv_vcf               // channel: [ val(meta), path(vcf) ]
        candidate_sv_vcf_tbi           = MANTA.out.candidate_sv_vcf_tbi           // channel: [ val(meta), path(tbi) ]
        diploid_sv_vcf                 = MANTA.out.diploid_sv_vcf                 // channel: [ val(meta), path(vcf) ]
        diploid_sv_vcf_tbi             = MANTA.out.diploid_sv_vcf_tbi             // channel: [ val(meta), path(tbi) ]
        filtered_diploid_sv_vcf        = BCFTOOLS_VIEW_MANTA.out.vcf              // channel: [ val(meta), path(vcf) ]
        versions                       = ch_versions
}
