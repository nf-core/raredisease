//
// Map to reference, fetch stats for each demultiplexed read pair, merge, mark duplicates, and index.
//

params.bwamem2_mem_options = [:]
params.samtools_idx_options = [:]
params.samtools_sort_options = [:]
params.samtools_stats_options = [:]
params.samtools_merge_options = [:]
params.markduplicates_options = [:]
params.samtools_idx_md_options = [:]

include { BWAMEM2_MEM } from '../../modules/nf-core/modules/bwamem2/mem/main'  addParams( options: params.bwamem2_mem_options )
include { SAMTOOLS_INDEX } from '../../modules/nf-core/modules/samtools/index/main' addParams( options: params.samtools_idx_options )
include { SAMTOOLS_SORT } from '../../modules/nf-core/modules/samtools/sort/main' addParams(options: params.samtools_sort_options )
include { SAMTOOLS_STATS } from '../../modules/nf-core/modules/samtools/stats/main' addParams(options: params.samtools_stats_options )
include { SAMTOOLS_MERGE } from '../../modules/nf-core/modules/samtools/merge/main' addParams( options: params.samtools_merge_options )
include { PICARD_MARKDUPLICATES as MARKDUPLICATES } from '../../modules/nf-core/modules/picard/markduplicates/main' addParams(options: params.markduplicates_options )


workflow MAPPING {
    take:
        reads_input // channel: [ val(meta), reads_input ]
        index       // channel: /path/to/bwamem2/index/

    main:
        // Map, sort, and index
        BWAMEM2_MEM ( reads_input, index )
        SAMTOOLS_SORT ( BWAMEM2_MEM.out.bam )
        SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )

        // Join the mapped bam + bai paths by their keys for stats
        // Get stats for each demultiplexed read pair.
        bam_sorted_indexed = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)
        SAMTOOLS_STATS ( bam_sorted_indexed )

        // Merge multiple lane samples and index
        SAMTOOLS_SORT.out.bam
        .map{ meta, bam ->
            new_meta = meta.clone()                 // clone to avoid overriding the global meta
            new_meta.id = new_meta.id.split('_')[0] // access the .id attribute of meta to split samplename_lane into samplename
            [new_meta, bam]}                        // end the closure to return newly modified channel
        .groupTuple(by: 0)                          // group them bam paths with the same [ [samplename], [bam path, bam path, ..] ]
        .branch{                                    // branch the channel into multiple channels (single, multiple) depending on size of list
            single: it[1].size() == 1
            multiple: it[1].size() > 1
            }
        .set{ bams }                                // create a new multi-channel named bams

        bams.multiple.view()

        // TODO: If there are no samples to merge, skip the process
        SAMTOOLS_MERGE ( bams.multiple )
        prepared_bam = bams.single.mix(SAMTOOLS_MERGE.out.bam)

        // Marking duplicates
        MARKDUPLICATES ( prepared_bam )

    emit:
        stats                  = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), [ stats ] ]
        metrics                = MARKDUPLICATES.out.metrics     // channel: [ val(meta), [ metrics ] ]
        marked_bam             = MARKDUPLICATES.out.bam         // channel: [ val(meta), [ marked_bam ] ]
        marked_bai             = MARKDUPLICATES.out.bai         // channel: [ val(meta), [ marked_bai ] ]

        // Collect versions
        bwamem2_version        = BWAMEM2_MEM.out.version        //      path: *.version.txt
        picard_version         = MARKDUPLICATES.out.version     //      path: *.version.txt
        samtools_version       = SAMTOOLS_SORT.out.version      //      path: *.version.txt
}
