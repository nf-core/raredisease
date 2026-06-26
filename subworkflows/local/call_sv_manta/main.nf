//
// A structural variant caller workflow for manta
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_MANTA } from '../../../modules/nf-core/bcftools/view/main.nf'
include { MANTA_GERMLINE as MANTA              } from '../../../modules/nf-core/manta/germline/main'

workflow CALL_SV_MANTA {
    take:
        ch_bam          // channel: [mandatory] [ val(meta), path(bam) ]
        ch_bai          // channel: [mandatory] [ val(meta), path(bai) ]
        ch_genome_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai   // channel: [mandatory] [ val(meta), path(fai) ]
        ch_case_info    // channel: [mandatory] [ val(case_info) ]
        ch_regions      // channel: [mandatory] [ path(bed), path(tbi) ]

    main:
        ch_bam.map{ _meta, bam -> bam }
            .collect(sort: { a, b -> a.getName() <=> b.getName() })
            .toList()
            .set { bam_file_list }

        ch_bai.map{ _meta, bai -> bai }
            .collect(sort: { a, b -> a.getName() <=> b.getName() })
            .toList()
            .set { bai_file_list }

        ch_case_info.combine(bam_file_list)
            .combine(bai_file_list)
            .combine(ch_regions)
            .set { manta_input }
        MANTA ( manta_input, ch_genome_fasta, ch_genome_fai, [] )

        MANTA.out.diploid_sv_vcf
            .join(MANTA.out.diploid_sv_vcf_tbi)
            .set {ch_filter_in}
        BCFTOOLS_VIEW_MANTA (ch_filter_in, [], [], [])

    emit:
        candidate_small_indels_vcf     = MANTA.out.candidate_small_indels_vcf     // channel: [ val(meta), path(vcf) ]
        candidate_small_indels_vcf_tbi = MANTA.out.candidate_small_indels_vcf_tbi // channel: [ val(meta), path(tbi) ]
        candidate_sv_vcf               = MANTA.out.candidate_sv_vcf               // channel: [ val(meta), path(vcf) ]
        candidate_sv_vcf_tbi           = MANTA.out.candidate_sv_vcf_tbi           // channel: [ val(meta), path(tbi) ]
        diploid_sv_vcf                 = MANTA.out.diploid_sv_vcf                 // channel: [ val(meta), path(vcf) ]
        diploid_sv_vcf_tbi             = MANTA.out.diploid_sv_vcf_tbi             // channel: [ val(meta), path(tbi) ]
        filtered_diploid_sv_vcf        = BCFTOOLS_VIEW_MANTA.out.vcf              // channel: [ val(meta), path(vcf) ]
}
