//
// A subworkflow to subsample MT alignments
//

include { BEDTOOLS_GENOMECOV      } from '../../modules/nf-core/bedtools/genomecov/main'
include { CALCULATE_SEED_FRACTION } from '../../modules/local/calculate_seed_fraction'
include { SAMTOOLS_VIEW           } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX          } from '../../modules/nf-core/samtools/index/main'

workflow SUBSAMPLE_MT {

    take:
        ch_mt_bam_bai          // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        val_mt_subsample_rd    // channel: [mandatory] [ val(read_dept) ]
        val_mt_subsample_seed  // channel: [mandatory] [ val(seed) ]

    main:
        ch_versions = Channel.empty()

        ch_mt_bam_bai.map {meta, bam, bai -> return [meta, bam, -1]}.set {ch_genomecov_in}

        BEDTOOLS_GENOMECOV (ch_genomecov_in, [], "genomecov")

        CALCULATE_SEED_FRACTION (
            BEDTOOLS_GENOMECOV.out.genomecov,
            val_mt_subsample_rd,
            val_mt_subsample_seed
        )
        .csv
        .join(ch_mt_bam_bai, failOnMismatch:true)
        .map{meta, seedfrac, bam, bai ->
            return [meta + [seedfrac: file(seedfrac).text.readLines()[0]], bam, bai]
        }
        .set { ch_subsample_in }

        SAMTOOLS_VIEW(ch_subsample_in, [[:],[]], [])

        SAMTOOLS_INDEX(SAMTOOLS_VIEW.out.bam)

        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())
        ch_versions = ch_versions.mix(CALCULATE_SEED_FRACTION.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
        versions = ch_versions  // channel: [ path(versions.yml) ]
}
