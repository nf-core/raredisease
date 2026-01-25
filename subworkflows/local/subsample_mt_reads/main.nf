//
// A subworkflow to subsample MT alignments
//

include { SAMTOOLS_VIEW               } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as SAM_TO_BAM } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_COLLATE            } from '../../../modules/nf-core/samtools/collate/main'
include { SAMTOOLS_SORT               } from '../../../modules/nf-core/samtools/sort/main'
include { GAWK                        } from '../../../modules/nf-core/gawk/main'

workflow SUBSAMPLE_MT_READS {

    take:
        ch_mt_bam_bai             // channel: [mandatory] [ val(meta), path(bam), path(bai) ]

    main:
        ch_versions = channel.empty()

        SAMTOOLS_VIEW(ch_mt_bam_bai, [[:],[]], [], '')

        SAMTOOLS_COLLATE(SAMTOOLS_VIEW.out.bam, [[:],[]])

        GAWK(SAMTOOLS_COLLATE.out.sam, [], false)

        GAWK.out.output.map {meta, sam -> return [meta, sam, []] }.set {ch_convert_to_bam}

        SAM_TO_BAM(ch_convert_to_bam, [[:],[]], [], '')

        SAMTOOLS_SORT(SAM_TO_BAM.out.bam, [[:],[]], 'bai')

        ch_versions = ch_versions.mix(SAMTOOLS_COLLATE.out.versions)
        ch_versions = ch_versions.mix(GAWK.out.versions)

    emit:
        versions = ch_versions  // channel: [ path(versions.yml) ]
}
