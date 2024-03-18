//
// A preprocessing workflow for Gens
//

include { GATK4_COLLECTREADCOUNTS as COLLECTREADCOUNTS        } from '../../modules/nf-core/gatk4/collectreadcounts/main'
include { GATK4_DENOISEREADCOUNTS as DENOISEREADCOUNTS_FEMALE } from '../../modules/nf-core/gatk4/denoisereadcounts/main'
include { GATK4_DENOISEREADCOUNTS as DENOISEREADCOUNTS_MALE   } from '../../modules/nf-core/gatk4/denoisereadcounts/main'
include { GENS as GENS_GENERATE                               } from '../../modules/local/gens/main'

workflow GENS {
    take:
        ch_bam_bai           // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_gvcf              // channel: [mandatory] [ val(meta), path(gvcf) ]
        ch_genome_fasta      // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai        // channel: [mandatory] [ val(meta), path(fai) ]
        ch_interval_list     // channel: [mandatory] [ path(interval_list) ]
        ch_pon_female        // channel: [mandatory] [ path(pon) ]
        ch_pon_male          // channel: [mandatory] [ path(pon) ]
        ch_gnomad_pos        // channel: [mandatory] [ path(gnomad_pos) ]
        ch_case_info         // channel: [mandatory] [ val(case_info) ]
        ch_genome_dictionary // channel: [mandatory] [ val(meta), path(dict) ]

    main:
        ch_versions = Channel.empty()

        ch_bam_bai
            .combine(ch_interval_list)
            .set { ch_bam_bai_intervals }

        COLLECTREADCOUNTS (
            ch_bam_bai_intervals,
            ch_genome_fasta,
            ch_genome_fai,
            ch_genome_dictionary
        )

        COLLECTREADCOUNTS.out.hdf5
            .branch { meta, counts ->
                female: meta.sex.equals(2) || meta.sex.equals(0)
                male: meta.sex.equals(1)
            }
            .set { ch_denoisereadcounts_in }

        DENOISEREADCOUNTS_FEMALE (
            ch_denoisereadcounts_in.female,
            ch_pon_female
        )

        DENOISEREADCOUNTS_MALE (
            ch_denoisereadcounts_in.male,
            ch_pon_male
        )
        DENOISEREADCOUNTS_FEMALE.out.standardized
            .mix(DENOISEREADCOUNTS_MALE.out.standardized)
            .set { ch_denoisereadcounts_out }

        GENS_GENERATE (
            ch_denoisereadcounts_out,
            ch_gvcf,
            ch_gnomad_pos
        )

        ch_versions = ch_versions.mix(COLLECTREADCOUNTS.out.versions.first())
        ch_versions = ch_versions.mix(DENOISEREADCOUNTS_FEMALE.out.versions.first())
        ch_versions = ch_versions.mix(DENOISEREADCOUNTS_MALE.out.versions.first())
        ch_versions = ch_versions.mix(GENS_GENERATE.out.versions.first())

    emit:
        gens_cov_bed_gz = GENS_GENERATE.out.cov // channel: [ val(meta), path(bed) ]
        gens_baf_bed_gz = GENS_GENERATE.out.baf // channel: [ val(meta), path(bed) ]
        versions        = ch_versions           // channel: [ path(versions.yml) ]
}
