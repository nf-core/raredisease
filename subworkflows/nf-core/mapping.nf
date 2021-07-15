/*
========================================================================================
    MAPPING
========================================================================================
*/

params.bwamem2_idx_options = [:]
params.bwamem2_mem_options = [:]
params.samtools_idx_options = [:]
params.samtools_merge_options = [:]
params.samtools_sort_options = [:]
params.samtools_stats_options = [:]
params.samtools_view_options = [:]

include { BWAMEM2_INDEX } from '../../modules/nf-core/modules/bwamem2/index/main'  addParams( options: params.bwamem2_idx_options )
include { BWAMEM2_MEM } from '../../modules/nf-core/modules/bwamem2/mem/main'  addParams( options: params.bwamem2_mem_options )
include { SAMTOOLS_INDEX } from '../../modules/nf-core/modules/samtools/index/main' addParams(options: params.samtools_idx_options )
include { SAMTOOLS_MERGE } from '../../modules/nf-core/modules/samtools/merge/main' addParams(options: params.samtools_merge_options )
include { SAMTOOLS_SORT } from '../../modules/nf-core/modules/samtools/sort/main' addParams(options: params.samtools_sort_options )
include { SAMTOOLS_STATS } from '../../modules/nf-core/modules/samtools/stats/main' addParams(options: params.samtools_stats_options )
include { SAMTOOLS_VIEW } from '../../modules/nf-core/modules/samtools/view/main' addParams(options: params.samtools_view_options )


workflow MAPPING {
    take:
        reads_input // channel: [mandatory] meta, reads_input
        fasta // channel: [mandatory] fasta

    main:

    BWAMEM2_INDEX ( fasta )
    BWAMEM2_MEM ( reads_input, BWAMEM2_INDEX.out.index )
    bwamem2_version = BWAMEM2_MEM.out.version

    SAMTOOLS_SORT( BWAMEM2_MEM.out.bam )
    sorted_bam_bwa = SAMTOOLS_SORT.out.bam

    sorted_bam_bwa.map{ meta, bam ->
        meta.id = meta.id.split('_')[0]
        [meta, bam]
    }.groupTuple().branch{
        single: it[1].size() == 1
        multiple: it[1].size() > 1
    }.set{ bam_bwa }

    SAMTOOLS_MERGE( bam_bwa.multiple )
    SAMTOOLS_INDEX( SAMTOOLS_MERGE.out.merged_bam )
    samtools_version = SAMTOOLS_SORT.out.version

    emit:
    bam = SAMTOOLS_MERGE.out.merged_bam
    bai = SAMTOOLS_INDEX.out.bai

    bwamem2_version
    samtools_version
}

