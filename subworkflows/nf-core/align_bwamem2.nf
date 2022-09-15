//
// Map to reference, fetch stats for each demultiplexed read pair, merge, mark duplicates, and index.
//

include { BWAMEM2_MEM                             } from '../../modules/nf-core/modules/bwamem2/mem/main'
include { SAMTOOLS_INDEX                          } from '../../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MD     } from '../../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_STATS                          } from '../../modules/nf-core/modules/samtools/stats/main'
include { SAMTOOLS_MERGE                          } from '../../modules/nf-core/modules/samtools/merge/main'
include { PICARD_MARKDUPLICATES as MARKDUPLICATES } from '../../modules/nf-core/modules/picard/markduplicates/main'


workflow ALIGN_BWAMEM2 {
    take:
        reads_input // channel: [ val(meta), reads_input ]
        index       // channel: [ /path/to/bwamem2/index/ ]
        fasta       // channel: [genome.fasta]
        fai         // channel: [genome.fai]

    main:
        ch_versions = Channel.empty()

        // Map, sort, and index
        BWAMEM2_MEM ( reads_input, index, true )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

        SAMTOOLS_INDEX ( BWAMEM2_MEM.out.bam )

        // Join the mapped bam + bai paths by their keys for stats
        // Get stats for each demultiplexed read pair.
        bam_sorted_indexed = BWAMEM2_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai)
        SAMTOOLS_STATS ( bam_sorted_indexed, [] )

        // Merge multiple lane samples and index
        BWAMEM2_MEM.out.bam
        .map{ meta, bam ->
            new_meta = meta.clone()                 // clone to avoid overriding the global meta
            new_meta.id = new_meta.id.split('_')[0] // access the .id attribute of meta to split samplename_lane into samplename
            new_meta.read_group = []                // remove read group from cloned meta to get a single file per sample
            [new_meta, bam]}                        // end the closure to return newly modified channel
        .groupTuple(by: 0)                          // group them bam paths with the same [ [samplename], [bam path, bam path, ..] ]
        .branch{                                    // branch the channel into multiple channels (single, multiple) depending on size of list
            single: it[1].size() == 1
            multiple: it[1].size() > 1
            }
        .set{ bams }                                // create a new multi-channel named bams

        // TODO: If there are no samples to merge, skip the process
        SAMTOOLS_MERGE ( bams.multiple, fasta, fai )
        prepared_bam = bams.single.mix(SAMTOOLS_MERGE.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

        // Marking duplicates
        MARKDUPLICATES ( prepared_bam )
        SAMTOOLS_INDEX_MD ( MARKDUPLICATES.out.bam )
        ch_versions = ch_versions.mix(MARKDUPLICATES.out.versions)

    emit:
        stats                  = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), [ stats ] ]
        metrics                = MARKDUPLICATES.out.metrics     // channel: [ val(meta), [ metrics ] ]
        marked_bam             = MARKDUPLICATES.out.bam         // channel: [ val(meta), [ marked_bam ] ]
        marked_bai             = SAMTOOLS_INDEX_MD.out.bai      // channel: [ val(meta), [ marked_bai ] ]

        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
