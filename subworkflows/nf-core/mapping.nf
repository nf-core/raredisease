//
// Map to reference, sort + index + stats each alignment file.
//

params.bwamem2_idx_options = [:]
params.bwamem2_mem_options = [:]
params.samtools_idx_options = [:]
params.samtools_sort_options = [:]
params.samtools_stats_options = [:]

include { BWAMEM2_INDEX } from '../../modules/nf-core/modules/bwamem2/index/main'  addParams( options: params.bwamem2_idx_options )
include { BWAMEM2_MEM } from '../../modules/nf-core/modules/bwamem2/mem/main'  addParams( options: params.bwamem2_mem_options )
include { SAMTOOLS_INDEX } from '../../modules/nf-core/modules/samtools/index/main' addParams(options: params.samtools_idx_options )
include { SAMTOOLS_SORT } from '../../modules/nf-core/modules/samtools/sort/main' addParams(options: params.samtools_sort_options )
include { SAMTOOLS_STATS } from '../../modules/nf-core/modules/samtools/stats/main' addParams(options: params.samtools_stats_options )


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

        // Join the mapped bam + bai paths by their keys
        bam_sorted_indexed = Channel.empty()
        bam_sorted_indexed = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)

        // Get stats
        SAMTOOLS_STATS ( bam_sorted_indexed )
        // TODO: MIP adds an add'l line to the stats file: percentage mapped: ((mapped/total raw) * 100). I think we should write a process to handle this.

        // Collect versions
        bwamem2_version = BWAMEM2_MEM.out.version
        samtools_version = SAMTOOLS_SORT.out.version

    emit:
        bam = SAMTOOLS_SORT.out.bam
        bai = SAMTOOLS_INDEX.out.bai
        stats = SAMTOOLS_STATS.out.stats

        bwamem2_version
        samtools_version
}
