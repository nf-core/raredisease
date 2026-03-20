//
// A subworkflow to annotate structural variants.
//

include { SENTIEON_BWAMEM                          } from '../../../modules/nf-core/sentieon/bwamem/main'
include { SENTIEON_DATAMETRICS                     } from '../../../modules/nf-core/sentieon/datametrics/main'
include { SENTIEON_DEDUP                           } from '../../../modules/nf-core/sentieon/dedup/main'
include { SENTIEON_READWRITER                      } from '../../../modules/nf-core/sentieon/readwriter/main'
include { SAMTOOLS_VIEW as EXTRACT_ALIGNMENTS      } from '../../../modules/nf-core/samtools/view/main'

workflow ALIGN_SENTIEON {
    take:
        ch_bwa_index           // channel: [mandatory] [ val(meta), path(bwa_index) ]
        ch_genome_fai          // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_input_reads         // channel: [mandatory] [ val(meta), path(reads_input) ]
        val_extract_alignments //  string: boolean
        val_platform           //  string: [mandatory] default: illumina

    main:

        SENTIEON_BWAMEM ( ch_input_reads, ch_bwa_index, ch_genome_fasta, ch_genome_fai )

        SENTIEON_BWAMEM.out
            .bam_and_bai
            .map{ meta, bam, bai ->
                def new_id   = meta.sample
                def new_meta = meta + [id:new_id, read_group:"\'@RG\\tID:" + new_id + "\\tPL:" + val_platform + "\\tSM:" + new_id + "\'"] - meta.subMap('lane')
                [groupKey(new_meta, new_meta.num_lanes), bam, bai]
                }
            .groupTuple()
            .branch{ _meta, bam, _bai ->
                single: bam.size() == 1
                multiple: bam.size() > 1
                }
            .set{ merge_bams_in }

        SENTIEON_READWRITER ( merge_bams_in.multiple, ch_genome_fasta, ch_genome_fai )
        ch_bam_bai = merge_bams_in.single.mix(SENTIEON_READWRITER.out.output_index)

        // GET ALIGNMENT FROM SELECTED CONTIGS
        if (val_extract_alignments) {
            EXTRACT_ALIGNMENTS( ch_bam_bai, ch_genome_fasta.join(ch_genome_fai).collect(), [], 'bai')
            ch_bam_bai = EXTRACT_ALIGNMENTS.out.bam.join(EXTRACT_ALIGNMENTS.out.bai, failOnMismatch:true, failOnDuplicate:true)
        }

        SENTIEON_DATAMETRICS ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, false )

        SENTIEON_DEDUP ( ch_bam_bai, ch_genome_fasta, ch_genome_fai )

        ch_publish = SENTIEON_DEDUP.out.bam
            .mix(SENTIEON_DEDUP.out.bai)
            .mix(SENTIEON_DEDUP.out.score)
            .mix(SENTIEON_DEDUP.out.metrics)
            .mix(SENTIEON_DEDUP.out.metrics_multiqc_tsv)
            .map { meta, value -> ['alignment/', [meta, value]] }

    emit:
        marked_bam  = SENTIEON_DEDUP.out.bam                             // channel: [ val(meta), path(bam) ]
        marked_bai  = SENTIEON_DEDUP.out.bai                             // channel: [ val(meta), path(bai) ]
        mq_metrics  = SENTIEON_DATAMETRICS.out.mq_metrics.ifEmpty(null)  // channel: [ val(meta), path(mq_metrics) ]
        qd_metrics  = SENTIEON_DATAMETRICS.out.qd_metrics.ifEmpty(null)  // channel: [ val(meta), path(qd_metrics) ]
        gc_metrics  = SENTIEON_DATAMETRICS.out.gc_metrics.ifEmpty(null)  // channel: [ val(meta), path(gc_metrics) ]
        gc_summary  = SENTIEON_DATAMETRICS.out.gc_summary.ifEmpty(null)  // channel: [ val(meta), path(gc_summary) ]
        aln_metrics = SENTIEON_DATAMETRICS.out.aln_metrics.ifEmpty(null) // channel: [ val(meta), path(aln_metrics) ]
        is_metrics  = SENTIEON_DATAMETRICS.out.is_metrics.ifEmpty(null)  // channel: [ val(meta), path(is_metrics) ]
        publish = ch_publish                                              // channel: [ val(destination), val(value) ]
}
