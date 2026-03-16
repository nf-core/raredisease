//
// A preprocessing workflow for Gens
//

include { GATK4_COLLECTREADCOUNTS as COLLECTREADCOUNTS        } from '../../../modules/nf-core/gatk4/collectreadcounts/main'
include { GATK4_DENOISEREADCOUNTS as DENOISEREADCOUNTS_FEMALE } from '../../../modules/nf-core/gatk4/denoisereadcounts/main'
include { GATK4_DENOISEREADCOUNTS as DENOISEREADCOUNTS_MALE   } from '../../../modules/nf-core/gatk4/denoisereadcounts/main'
include { PREPARECOVANDBAF as GENS_GENERATE                   } from '../../../modules/nf-core/gens/preparecovandbaf/main'

workflow GENS {
    take:
        ch_bam_bai           // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_dictionary // channel: [mandatory] [ val(meta), path(dict) ]
        ch_genome_fai        // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_fasta      // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_gnomad_pos        // channel: [mandatory] [ path(gnomad_pos) ]
        ch_gvcf              // channel: [mandatory] [ val(meta), path(gvcf) ]
        ch_gvcf_tbi          // channel: [mandatory] [ val(meta), path(gvcf.tbi) ]
        ch_interval_list     // channel: [mandatory] [ path(interval_list) ]
        ch_pon_female        // channel: [mandatory] [ path(pon) ]
        ch_pon_male          // channel: [mandatory] [ path(pon) ]

    main:
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
            .branch { meta, _counts ->
                female: meta.sex.toString().matches('2|other|0')
                male: meta.sex == 1
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
            ch_denoisereadcounts_out
                .join(ch_gvcf)
                .join(ch_gvcf_tbi),
            ch_gnomad_pos
        )

    emit:
        gens_baf_bed_gz  = GENS_GENERATE.out.baf_gz  // channel: [ val(meta), path(bed.gz) ]
        gens_baf_bed_tbi = GENS_GENERATE.out.baf_tbi // channel: [ val(meta), path(tbi) ]
        gens_cov_bed_gz  = GENS_GENERATE.out.cov_gz  // channel: [ val(meta), path(bed.gz) ]
        gens_cov_bed_tbi = GENS_GENERATE.out.cov_tbi // channel: [ val(meta), path(tbi) ]
}
