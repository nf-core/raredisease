//
// Map to reference, fetch stats for each demultiplexed read pair, merge, mark duplicates, and index.
//

include { BWAMEM2_MEM                              } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWAMEME_MEM                              } from '../../../modules/nf-core/bwameme/mem/main'
include { BWA_MEM as BWA                           } from '../../../modules/nf-core/bwa/mem/main'
include { FASTDUP                                  } from '../../../modules/nf-core/fastdup/main'
include { PICARD_MARKDUPLICATES as MARKDUPLICATES  } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_ALIGN   } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_EXTRACT } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MARKDUP } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                           } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_STATS                           } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_VIEW as EXTRACT_ALIGNMENTS      } from '../../../modules/nf-core/samtools/view/main'


workflow ALIGN_BWA_BWAMEM2_BWAMEME {
    take:
        ch_bwa_index           // channel: [mandatory] [ val(meta), path(bwa_index) ]
        ch_bwamem2_index       // channel: [mandatory] [ val(meta), path(bwamem2_index) ]
        ch_bwameme_index       // channel: [mandatory] [ val(meta), path(bwameme_index) ]
        ch_genome_fai          // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_input_reads         // channel: [mandatory] [ val(meta), path(reads_input) ]
        val_aligner            // string:  'bwa', 'bwamem2', 'bwameme', or 'sentieon'
        val_duplicates_marker  // string:  'markduplicates' or 'fastdup', default: 'markduplicates'
        val_extract_alignments // boolean
        val_mbuffer_mem        // integer: [mandatory] default: 3072
        val_platform           // string:  [mandatory] default: illumina
        val_sort_threads       // integer: [mandatory] default: 4

    main:
        // Map, sort, and index
        if (val_aligner.equals("bwa")) {
            BWA ( ch_input_reads, ch_bwa_index, ch_genome_fasta, true )
            ch_align = BWA.out.bam
        } else if (val_aligner.equals("bwameme")) {
            BWAMEME_MEM ( ch_input_reads, ch_bwameme_index, ch_genome_fasta, true, val_mbuffer_mem, val_sort_threads )
            ch_align = BWAMEME_MEM.out.bam
        } else {
            BWAMEM2_MEM ( ch_input_reads, ch_bwamem2_index, ch_genome_fasta, true )
            ch_align    = BWAMEM2_MEM.out.bam
        }

        SAMTOOLS_INDEX_ALIGN ( ch_align )

        // Get stats for each demultiplexed read pair.
        bam_sorted_indexed = ch_align.join(SAMTOOLS_INDEX_ALIGN.out.bai, failOnMismatch:true, failOnDuplicate:true)
        SAMTOOLS_STATS ( bam_sorted_indexed, [[],[]] )

        // Merge multiple lane samples and index
        ch_align
            .map{ meta, bam ->
                    def new_id   = meta.sample
                    def new_meta = meta + [id:new_id, read_group:"\'@RG\\tID:" + new_id + "\\tPL:" + val_platform + "\\tSM:" + new_id + "\'"] - meta.subMap('lane','data_type')
                    [groupKey(new_meta, new_meta.num_lanes), bam]
                }
            .groupTuple()
            .branch{ _meta, bam ->
                single: bam.size() == 1
                multiple: bam.size() > 1
                }
            .set{ bams }

        // If there are no samples to merge, skip the process
        SAMTOOLS_MERGE ( bams.multiple.map { it -> it + [[]] }, ch_genome_fasta.join(ch_genome_fai).map{meta,fasta,fai-> return [meta,fasta,fai,[]]}.collect())
        prepared_bam = bams.single.mix(SAMTOOLS_MERGE.out.bam)

        // GET ALIGNMENT FROM SELECTED CONTIGS
        if (val_extract_alignments) {
            SAMTOOLS_INDEX_EXTRACT ( prepared_bam )
            extract_bam_sorted_indexed = prepared_bam.join(SAMTOOLS_INDEX_EXTRACT.out.bai, failOnMismatch:true, failOnDuplicate:true)
            EXTRACT_ALIGNMENTS( extract_bam_sorted_indexed, ch_genome_fasta.join(ch_genome_fai).collect(), [], '')
            prepared_bam = EXTRACT_ALIGNMENTS.out.bam
        }

        // Marking duplicates
        if (val_duplicates_marker == "markduplicates") {
            MARKDUPLICATES ( prepared_bam, ch_genome_fasta, ch_genome_fai )
            SAMTOOLS_INDEX_MARKDUP (MARKDUPLICATES.out.bam)

            ch_marked_bam = MARKDUPLICATES.out.bam
            ch_marked_bai = SAMTOOLS_INDEX_MARKDUP.out.bai
            ch_marked_csi = SAMTOOLS_INDEX_MARKDUP.out.csi
            ch_metrics    = MARKDUPLICATES.out.metrics
        } else {
            FASTDUP ( prepared_bam ) 
            ch_marked_bam = FASTDUP.out.bam
            ch_marked_bai = FASTDUP.out.bai
            ch_marked_csi = FASTDUP.out.csi
            ch_metrics    = FASTDUP.out.metrics
        }

        ch_publish = ch_marked_bam
            .mix(ch_metrics)
            .mix(ch_marked_bai)
            .mix(ch_marked_csi)
            .map {meta, value -> ['alignment/', [meta, value]] }  


    emit:

        marked_bam = ch_marked_bam // = MARKDUPLICATES.out.bam or FASTDUP.out.bam           // channel: [ val(meta), path(bam) ]
        marked_bai = ch_marked_bai  // = SAMTOOLS_INDEX_MARKDUP.out.bai or FASTDUP.out.bai   // channel: [ val(meta), path(bai) ]
        marked_csi = ch_marked_csi  // = SAMTOOLS_INDEX_MARKDUP.out.csi or FASTDUP.out.csi   // channel: [ val(meta), path(csi) ]
        metrics    = ch_metrics     // = MARKDUPLICATES.out.metrics or FASTDUP.out.metrics   // channel: [ val(meta), path(metrics) ]
        stats      = SAMTOOLS_STATS.out.stats                              // channel: [ val(meta), path(stats) ]
        publish    = ch_publish                                                 // channel: [ val(destination), val(value) ]
}
