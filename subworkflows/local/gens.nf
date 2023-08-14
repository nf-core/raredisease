//
// A preprocessing workflow for Gens
//

include { GATK4_COLLECTREADCOUNTS as COLLECTREADCOUNTS } from '../../modules/nf-core/gatk4/collectreadcounts/main'
include { GATK4_DENOISEREADCOUNTS as DENOISEREADCOUNTS } from '../../modules/nf-core/gatk4/denoisereadcounts/main'
include { GENS as GENS_GENERATE                        } from '../../modules/local/gens/main'

workflow GENS {
    take:
        ch_bam_bai            // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_vcf                // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_genome_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai         // channel: [mandatory] [ val(meta), path(fai) ]
        ch_interval_list      // channel: [mandatory] [ path(interval_list) ]
        ch_pon                // channel: [mandatory] [ path(pon) ]
        ch_gnomad_pos         // channel: [mandatory] [ path(gnomad_pos) ]
        ch_case_info          // channel: [mandatory] [ val(case_info) ]
        ch_genome_dictionary  // channel: [mandatory] [ val(meta), path(dict) ]

    main:
        ch_versions = Channel.empty()

        COLLECTREADCOUNTS (ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_sequence_dictionary, ch_interval_list)

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
