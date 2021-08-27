//
// Map to reference, fetch stats for each demultiplexed read pair, merge, mark duplicates, and index.
//

// params.bwamem2_idx_options = [:]
params.bwamem2_mem_options = [:]
params.samtools_idx_options = [:]
params.samtools_sort_options = [:]
params.samtools_stats_options = [:]
params.samtools_merge_options = [:]
params.markduplicates_options = [:]
params.samtools_idx_md_options = [:]

// include { BWAMEM2_INDEX } from '../../modules/nf-core/modules/bwamem2/index/main'  addParams( options: params.bwamem2_idx_options )
include { BWAMEM2_MEM } from '../../modules/nf-core/modules/bwamem2/mem/main'  addParams( options: params.bwamem2_mem_options )
include { SAMTOOLS_INDEX } from '../../modules/nf-core/modules/samtools/index/main' addParams(options: params.samtools_idx_options )
include { SAMTOOLS_SORT } from '../../modules/nf-core/modules/samtools/sort/main' addParams(options: params.samtools_sort_options )
include { SAMTOOLS_STATS } from '../../modules/nf-core/modules/samtools/stats/main' addParams(options: params.samtools_stats_options )
include { SAMTOOLS_MERGE } from '../../modules/nf-core/modules/samtools/merge/main' addParams( options: params.samtools_merge_options )
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MD } from '../../modules/nf-core/modules/samtools/index/main' addParams(options: params.samtools_idx_md_options )
include { PICARD_MARKDUPLICATES as MARKDUPLICATES } from '../../modules/nf-core/modules/picard/markduplicates/main' addParams(options: params.markduplicates_options )


workflow MAPPING {
    take:
        reads_input // channel: [mandatory] meta, reads_input
        // fasta // channel: [mandatory] fasta
        index // channel: /path/to/bwamem2/index/

    main:
        // Map, sort, and index
        BWAMEM2_MEM ( reads_input, index )
        SAMTOOLS_SORT ( BWAMEM2_MEM.out.bam )
        SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )

        // Join the mapped bam + bai paths by their keys for stats
        // Get stats for each demultiplexed read pair.
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
        merged_bam = bams_to_merge.single.mix(SAMTOOLS_MERGE.out.bam)

        // Marking duplicates + index
        MARKDUPLICATES ( merged_bam )
        marked_bam = MARKDUPLICATES.out.bam

        SAMTOOLS_INDEX_MD ( marked_bam )
        marked_bai = SAMTOOLS_INDEX_MD.out.bai


        // Collect versions
        bwamem2_version = BWAMEM2_MEM.out.version
        markduplicates_version = MARKDUPLICATES.out.version
        samtools_version = SAMTOOLS_SORT.out.version

    emit:
        stats
        marked_bam
        marked_bai


        bwamem2_version
        markduplicates_version
        samtools_version
}
