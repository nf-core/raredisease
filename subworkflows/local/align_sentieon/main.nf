//
// A subworkflow to annotate structural variants.
//

include { SENTIEON_BWAMEM                          } from '../../../modules/nf-core/sentieon/bwamem/main'
include { SENTIEON_DATAMETRICS                     } from '../../../modules/nf-core/sentieon/datametrics/main'
include { SENTIEON_DEDUP                           } from '../../../modules/nf-core/sentieon/dedup/main'
include { SENTIEON_READWRITER                      } from '../../../modules/nf-core/sentieon/readwriter/main'
include { SAMTOOLS_VIEW as EXTRACT_ALIGNMENTS      } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_EXTRACT } from '../../../modules/nf-core/samtools/index/main'

workflow ALIGN_SENTIEON {
    take:
        ch_reads_input     // channel: [mandatory] [ val(meta), path(reads_input) ]
        ch_genome_fasta    // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai      // channel: [mandatory] [ val(meta), path(fai) ]
        ch_bwa_index       // channel: [mandatory] [ val(meta), path(bwa_index) ]
        val_platform       // string:  [mandatory] default: illumina

    main:
        ch_versions = Channel.empty()

        SENTIEON_BWAMEM ( ch_reads_input, ch_bwa_index, ch_genome_fasta, ch_genome_fai )

        SENTIEON_BWAMEM.out
            .bam_and_bai
            .map{ meta, bam, bai ->
                def new_id   = meta.sample
                def new_meta = meta + [id:new_id, read_group:"\'@RG\\tID:" + new_id + "\\tPL:" + val_platform + "\\tSM:" + new_id + "\'"] - meta.subMap('lane')
                [groupKey(new_meta, new_meta.num_lanes), bam, bai]
                }
            .groupTuple()
            .branch{
                single: it[1].size() == 1
                multiple: it[1].size() > 1
                }
            .set{ merge_bams_in }

        SENTIEON_READWRITER ( merge_bams_in.multiple, ch_genome_fasta, ch_genome_fai )
        ch_bam_bai = merge_bams_in.single.mix(SENTIEON_READWRITER.out.output_index)

        // GET ALIGNMENT FROM SELECTED CONTIGS
        if (params.extract_alignments) {
            EXTRACT_ALIGNMENTS( ch_bam_bai, ch_genome_fasta, [])
            ch_bam_bai = EXTRACT_ALIGNMENTS.out.bam
            SAMTOOLS_INDEX_EXTRACT ( EXTRACT_ALIGNMENTS.out.bam )
            ch_bam_bai = EXTRACT_ALIGNMENTS.out.bam.join(SAMTOOLS_INDEX_EXTRACT.out.bai, failOnMismatch:true, failOnDuplicate:true)
            ch_versions = ch_versions.mix(EXTRACT_ALIGNMENTS.out.versions.first())
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX_EXTRACT.out.versions.first())

        }

        SENTIEON_DATAMETRICS ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, false )

        SENTIEON_DEDUP ( ch_bam_bai, ch_genome_fasta, ch_genome_fai )

        ch_versions = ch_versions.mix(SENTIEON_BWAMEM.out.versions.first())
        ch_versions = ch_versions.mix(SENTIEON_DATAMETRICS.out.versions.first())
        ch_versions = ch_versions.mix(SENTIEON_READWRITER.out.versions.first())
        ch_versions = ch_versions.mix(SENTIEON_DEDUP.out.versions.first())

    emit:
        marked_bam  = SENTIEON_DEDUP.out.bam                             // channel: [ val(meta), path(bam) ]
        marked_bai  = SENTIEON_DEDUP.out.bai                             // channel: [ val(meta), path(bai) ]
        mq_metrics  = SENTIEON_DATAMETRICS.out.mq_metrics.ifEmpty(null)  // channel: [ val(meta), path(mq_metrics) ]
        qd_metrics  = SENTIEON_DATAMETRICS.out.qd_metrics.ifEmpty(null)  // channel: [ val(meta), path(qd_metrics) ]
        gc_metrics  = SENTIEON_DATAMETRICS.out.gc_metrics.ifEmpty(null)  // channel: [ val(meta), path(gc_metrics) ]
        gc_summary  = SENTIEON_DATAMETRICS.out.gc_summary.ifEmpty(null)  // channel: [ val(meta), path(gc_summary) ]
        aln_metrics = SENTIEON_DATAMETRICS.out.aln_metrics.ifEmpty(null) // channel: [ val(meta), path(aln_metrics) ]
        is_metrics  = SENTIEON_DATAMETRICS.out.is_metrics.ifEmpty(null)  // channel: [ val(meta), path(is_metrics) ]
        versions    = ch_versions                                        // channel: [ path(versions.yml) ]
}
