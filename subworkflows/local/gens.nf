//
// A preprocessing workflow for Gens
//

include { GATK4_COLLECTREADCOUNTS as COLLECTREADCOUNTS } from '../../modules/local/gatk4/collectreadcounts/main'
include { GATK4_DENOISEREADCOUNTS as DENOISEREADCOUNTS } from '../../modules/local/gatk4/denoisereadcounts/main'
include { GENS as GENS_GENERATE                        } from '../../modules/local/gens/main'

workflow GENS {
    take:
        ch_bam           // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_vcf           // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_fasta         // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_fai           // channel: [mandatory] [ path(fai) ]
        ch_interval_list // channel: [mandatory] [ path(interval_list) ]
        ch_pon           // channel: [mandatory] [ path(pon) ]
        ch_gnomad_pos    // channel: [mandatory] [ path(gnomad_pos) ]
        ch_case_info     // channel: [mandatory] [ val(case_info) ]
        ch_seq_dict      // channel: [mandatory] [ path(dict) ]

    main:
        ch_versions = Channel.empty()

        COLLECTREADCOUNTS (ch_bam, ch_fasta, ch_fai, ch_seq_dict, ch_interval_list)

        DENOISEREADCOUNTS (COLLECTREADCOUNTS.out.read_counts, ch_pon)

        GENS_GENERATE (DENOISEREADCOUNTS.out.standardized_read_counts, ch_vcf.map { meta, vcf -> vcf }, ch_gnomad_pos)

        ch_versions = ch_versions.mix(COLLECTREADCOUNTS.out.versions.first())
        ch_versions = ch_versions.mix(DENOISEREADCOUNTS.out.versions.first())
        ch_versions = ch_versions.mix(GENS_GENERATE.out.versions.first())

    emit:
        gens_cov_bed_gz = GENS_GENERATE.out.cov // channel: [ val(meta), path(bed) ]
        gens_baf_bed_gz = GENS_GENERATE.out.baf // channel: [ val(meta), path(bed) ]
        versions        = ch_versions           // channel: [ path(versions.yml) ]
}
