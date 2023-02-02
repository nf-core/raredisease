//
// Map to reference, fetch stats for each demultiplexed read pair, merge, mark duplicates, and index.
//

include { BWAMEM2_MEM                              } from '../../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALIGN   } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MARKDUP } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS                           } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_MERGE                           } from '../../../modules/nf-core/samtools/merge/main'
include { PICARD_MARKDUPLICATES as MARKDUPLICATES  } from '../../../modules/nf-core/picard/markduplicates/main'


workflow ALIGN_BWAMEM2 {
    take:
        reads_input // channel: [ val(meta), reads_input ]
        index       // channel: [ /path/to/bwamem2/index/ ]
        fasta       // channel: [genome.fasta]
        fai         // channel: [genome.fai]
        platform    // params.platform

    main:
        ch_versions = Channel.empty()

        // Map, sort, and index
        BWAMEM2_MEM ( reads_input, index, true )

        SAMTOOLS_INDEX_ALIGN ( BWAMEM2_MEM.out.bam )

        // Get stats for each demultiplexed read pair.
        bam_sorted_indexed = BWAMEM2_MEM.out.bam.join(SAMTOOLS_INDEX_ALIGN.out.bai)
        SAMTOOLS_STATS ( bam_sorted_indexed, [] )

        // Merge multiple lane samples and index
        BWAMEM2_MEM.out.bam
            .map{ meta, bam ->
                    new_meta            = meta.clone()
                    new_meta.id         = new_meta.id.split('_')[0]
                    new_meta.read_group = "\'@RG\\tID:" + new_meta.id + "\\tPL:" + platform + "\\tSM:" + new_meta.id + "\'"
                    [new_meta, bam]
                }
            .groupTuple(by: 0)
            .branch{
                single: it[1].size() == 1
                multiple: it[1].size() > 1
                }
            .set{ bams }

        // If there are no samples to merge, skip the process
        SAMTOOLS_MERGE ( bams.multiple, fasta, fai )
        prepared_bam = bams.single.mix(SAMTOOLS_MERGE.out.bam)

        // Marking duplicates
        MARKDUPLICATES ( prepared_bam , fasta, fai )
        SAMTOOLS_INDEX_MARKDUP ( MARKDUPLICATES.out.bam )

        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_ALIGN.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())
        ch_versions = ch_versions.mix(MARKDUPLICATES.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MARKDUP.out.versions.first())

    emit:
        stats                  = SAMTOOLS_STATS.out.stats
        metrics                = MARKDUPLICATES.out.metrics
        marked_bam             = MARKDUPLICATES.out.bam
        marked_bai             = SAMTOOLS_INDEX_MARKDUP.out.bai
        versions               = ch_versions.ifEmpty(null)
}
