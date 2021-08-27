//
// Map to reference, sort + index + stats each alignment file, and then merge.
//

params.bwamem2_idx_options = [:]
params.bwamem2_mem_options = [:]
params.samtools_idx_options = [:]
params.samtools_sort_options = [:]
params.samtools_stats_options = [:]
params.samtools_merge_options = [:]

include { BWAMEM2_INDEX } from '../../modules/nf-core/modules/bwamem2/index/main'  addParams( options: params.bwamem2_idx_options )
include { BWAMEM2_MEM } from '../../modules/nf-core/modules/bwamem2/mem/main'  addParams( options: params.bwamem2_mem_options )
include { SAMTOOLS_INDEX } from '../../modules/nf-core/modules/samtools/index/main' addParams(options: params.samtools_idx_options )
include { SAMTOOLS_SORT } from '../../modules/nf-core/modules/samtools/sort/main' addParams(options: params.samtools_sort_options )
include { SAMTOOLS_STATS } from '../../modules/nf-core/modules/samtools/stats/main' addParams(options: params.samtools_stats_options )
include { SAMTOOLS_MERGE } from '../../modules/nf-core/modules/samtools/merge/main' addParams( options: params.samtools_merge_options )
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_2 } from '../../modules/nf-core/modules/samtools/index/main' addParams(options: params.samtools_idx_options )


workflow MAPPING {
    take:
        reads_input // channel: [mandatory] meta, reads_input
        fasta // channel: [mandatory] fasta

    main:
        // Index
        BWAMEM2_INDEX ( fasta )

        // Map, sort, and index
        BWAMEM2_MEM ( reads_input, BWAMEM2_INDEX.out.index )
        SAMTOOLS_SORT ( BWAMEM2_MEM.out.bam )
        SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )

        // Join the mapped bam + bai paths by their keys for stats
        // Get stats for each demultiplexed read pair.
        bam_sorted_indexed = Channel.empty()
        bam_sorted_indexed = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)
        SAMTOOLS_STATS ( bam_sorted_indexed )
        stats = SAMTOOLS_STATS.out.stats

        // Merge multiple lane samples and index
        SAMTOOLS_SORT.out.bam.map{ meta, bam ->
            new_meta = meta.clone()
            new_meta.id = new_meta.id.split('_')[0]
            [new_meta, bam]
        }.groupTuple().branch{
            single: it[1].size() == 1
            multiple: it[1].size() > 1
        }.set{ bams_to_merge }

        SAMTOOLS_MERGE ( bams_to_merge.multiple )
        bam = bams_to_merge.single.mix(SAMTOOLS_MERGE.out.bam)

        SAMTOOLS_INDEX_2 ( bam )
        bai = SAMTOOLS_INDEX_2.out.bai


        // Collect versions
        bwamem2_version = BWAMEM2_MEM.out.version
        samtools_version = SAMTOOLS_SORT.out.version

    emit:
        stats
        bam
        bai

        bwamem2_version
        samtools_version
}
