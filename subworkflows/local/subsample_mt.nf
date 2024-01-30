//
// A subworkflow to subsample MT alignments
//

include { BEDTOOLS_GENOMECOV                    } from '../../modules/nf-core/bedtools/genomecov/main'

workflow SUBSAMPLE_MT {

    take:
        ch_mt_marked_bam                // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ]

    main:

        ch_mt_marked_bam.map {meta, bam -> return [meta, bam, []]}.set {ch_genomecov_in}

        BEDTOOLS_GENOMECOV (ch_genomecov_in, [], [])

        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

    emit:
        versions = ch_versions  // channel: [ path(versions.yml) ]
}
