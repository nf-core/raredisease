//
// Map to reference, fetch stats for each demultiplexed read pair, merge, mark duplicates, and index.
//

include { BWA_MEM                                  } from '../../../modules/nf-core/bwa/mem/main'
include { BWAMEM2_MEM                              } from '../../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALIGN   } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MARKDUP } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS                           } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_MERGE                           } from '../../../modules/nf-core/samtools/merge/main'
include { PICARD_MARKDUPLICATES as MARKDUPLICATES  } from '../../../modules/nf-core/picard/markduplicates/main'


workflow ALIGN_BWA_BWAMEM2 {
    take:
        ch_reads_input   // channel: [mandatory] [ val(meta), path(reads_input) ]
        ch_bwa_index     // channel: [mandatory] [ val(meta), path(bwamem2_index) ]
        ch_bwamem2_index // channel: [mandatory] [ val(meta), path(bwamem2_index) ]
        ch_genome_fasta  // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai    // channel: [mandatory] [ val(meta), path(fai) ]
        val_platform     // string:  [mandatory] default: illumina

    main:
        ch_versions = Channel.empty()

        // Map, sort, and index
        if (params.aligner.equals("bwa")) {
            BWA_MEM ( ch_reads_input, ch_bwa_index, true )
            ch_align = BWA_MEM.out.bam
            ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
        } else {
            BWAMEM2_MEM ( ch_reads_input, ch_bwamem2_index, true )
            ch_align = BWAMEM2_MEM.out.bam
            ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
        }

        SAMTOOLS_INDEX_ALIGN ( ch_align )

        // Get stats for each demultiplexed read pair.
        bam_sorted_indexed = ch_align.join(SAMTOOLS_INDEX_ALIGN.out.bai, failOnMismatch:true, failOnDuplicate:true)
        SAMTOOLS_STATS ( bam_sorted_indexed, [[],[]] )

        // Merge multiple lane samples and index
        ch_align
            .map{ meta, bam ->
                    new_id   = meta.sample
                    new_meta = meta + [id:new_id, read_group:"\'@RG\\tID:" + new_id + "\\tPL:" + val_platform + "\\tSM:" + new_id + "\'"]
                    [groupKey(new_meta, new_meta.num_lanes), bam]
                }
            .groupTuple()
            .branch{
                single: it[1].size() == 1
                multiple: it[1].size() > 1
                }
            .set{ bams }

        // If there are no samples to merge, skip the process
        SAMTOOLS_MERGE ( bams.multiple, ch_genome_fasta, ch_genome_fai )
        prepared_bam = bams.single.mix(SAMTOOLS_MERGE.out.bam)

        // Marking duplicates
        MARKDUPLICATES ( prepared_bam , ch_genome_fasta, ch_genome_fai )
        SAMTOOLS_INDEX_MARKDUP ( MARKDUPLICATES.out.bam )

        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_ALIGN.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())
        ch_versions = ch_versions.mix(MARKDUPLICATES.out.versions.first())
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_MARKDUP.out.versions.first())

    emit:
        stats       = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), path(stats) ]
        metrics     = MARKDUPLICATES.out.metrics     // channel: [ val(meta), path(metrics) ]
        marked_bam  = MARKDUPLICATES.out.bam         // channel: [ val(meta), path(bam) ]
        marked_bai  = SAMTOOLS_INDEX_MARKDUP.out.bai // channel: [ val(meta), path(bai) ]
        versions    = ch_versions                    // channel: [ path(versions.yml) ]
}
