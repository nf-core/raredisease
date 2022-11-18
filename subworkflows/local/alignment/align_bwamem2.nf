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
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

        SAMTOOLS_INDEX_ALIGN ( BWAMEM2_MEM.out.bam )

        // Join the mapped bam + bai paths by their keys for stats
        // Get stats for each demultiplexed read pair.
        bam_sorted_indexed = BWAMEM2_MEM.out.bam.join(SAMTOOLS_INDEX_ALIGN.out.bai)
        SAMTOOLS_STATS ( bam_sorted_indexed, [] )

        // Merge multiple lane samples and index
        BWAMEM2_MEM.out.bam
        .map{ meta, bam ->
            new_meta = meta.clone()                                                                                 // clone to avoid overriding the global meta
            new_meta.id = new_meta.id.split('_')[0]                                                                 // access the .id attribute of meta to split samplename_lane into samplename
            new_meta.read_group = "\'@RG\\tID:" + new_meta.id + "\\tPL:" + platform + "\\tSM:" + new_meta.id + "\'"
            [new_meta, bam]}                                                                                        // end the closure to return newly modified channel
        .groupTuple(by: 0)                                                                                          // group them bam paths with the same [ [samplename], [bam path, bam path, ..] ]
        .branch{                                                                                                    // branch the channel into multiple channels (single, multiple) depending on size of list
            single: it[1].size() == 1
            multiple: it[1].size() > 1
            }
        .set{ bams }                                // create a new multi-channel named bams

        bams.multiple.view()
        bams.single.view()
        // If there are no samples to merge, skip the process
        SAMTOOLS_MERGE ( bams.multiple, fasta, fai )
        prepared_bam = bams.single.mix(SAMTOOLS_MERGE.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

        // Marking duplicates
        MARKDUPLICATES ( prepared_bam , fasta, fai )
        SAMTOOLS_INDEX_MARKDUP ( MARKDUPLICATES.out.bam )
        ch_versions = ch_versions.mix(MARKDUPLICATES.out.versions)

    emit:
        stats                  = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), [ stats ] ]
        metrics                = MARKDUPLICATES.out.metrics     // channel: [ val(meta), [ metrics ] ]
        marked_bam             = MARKDUPLICATES.out.bam         // channel: [ val(meta), [ marked_bam ] ]
        marked_bai             = SAMTOOLS_INDEX_MARKDUP.out.bai      // channel: [ val(meta), [ marked_bai ] ]

        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
