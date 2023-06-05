//
// A subworkflow to annotate structural variants.
//

include { SENTIEON_BWAMEM         } from '../../../modules/local/sentieon/bwamem'
include { SENTIEON_DATAMETRICS    } from '../../../modules/local/sentieon/datametrics'
include { SENTIEON_LOCUSCOLLECTOR } from '../../../modules/local/sentieon/locuscollector'
include { SENTIEON_DEDUP          } from '../../../modules/local/sentieon/dedup'
include { SENTIEON_BQSR           } from '../../../modules/local/sentieon/bqsr'
include { SENTIEON_READWRITER     } from '../../../modules/local/sentieon/readwriter'

workflow ALIGN_SENTIEON {
    take:
        ch_reads_input     // channel: [mandatory] [ val(meta), path(reads_input) ]
        ch_genome_fasta    // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai      // channel: [mandatory] [ val(meta), path(fai) ]
        ch_bwa_index       // channel: [mandatory] [ val(meta), path(bwa_index) ]
        ch_known_dbsnp     // channel: [optional] [ path(known_dbsnp) ]
        ch_known_dbsnp_tbi // channel: [optional] [ path(known_dbsnp_tbi) ]
        val_platform       // string:  [mandatory] default: illumina

    main:
        ch_versions = Channel.empty()
        ch_bqsr_bam = Channel.empty()
        ch_bqsr_bai = Channel.empty()
        ch_bqsr_csv = Channel.empty()

        SENTIEON_BWAMEM ( ch_reads_input, ch_genome_fasta, ch_genome_fai, ch_bwa_index )

        SENTIEON_BWAMEM.out
            .bam
            .join(SENTIEON_BWAMEM.out.bai, failOnMismatch:true, failOnDuplicate:true)
            .map{ meta, bam, bai ->
                new_id   = meta.id.split('_')[0]
                new_meta = meta + [id:new_id, read_group:"\'@RG\\tID:" + new_id + "\\tPL:" + val_platform + "\\tSM:" + new_id + "\'"]
                [groupKey(new_meta, new_meta.num_lanes), bam, bai]
                }
            .groupTuple()
            .branch{
                single: it[1].size() == 1
                multiple: it[1].size() > 1
                }
            .set{ merge_bams_in }

        SENTIEON_READWRITER (merge_bams_in.multiple)
        ch_bam_bai = merge_bams_in.single.mix(SENTIEON_READWRITER.out.bam_bai)

        SENTIEON_DATAMETRICS (ch_bam_bai, ch_genome_fasta, ch_genome_fai )

        SENTIEON_LOCUSCOLLECTOR ( ch_bam_bai )

        ch_bam_bai
            .join(SENTIEON_LOCUSCOLLECTOR.out.score, failOnMismatch:true, failOnDuplicate:true)
            .join(SENTIEON_LOCUSCOLLECTOR.out.score_idx, failOnMismatch:true, failOnDuplicate:true)
            .set { ch_bam_bai_score }

        SENTIEON_DEDUP ( ch_bam_bai_score, ch_genome_fasta, ch_genome_fai )

        if (params.variant_caller == "sentieon") {
            SENTIEON_DEDUP.out.bam
                .join(SENTIEON_DEDUP.out.bai, failOnMismatch:true, failOnDuplicate:true)
                .set { ch_dedup_bam_bai }
            SENTIEON_BQSR ( ch_dedup_bam_bai, ch_genome_fasta, ch_genome_fai, ch_known_dbsnp, ch_known_dbsnp_tbi )
            ch_bqsr_bam = SENTIEON_BQSR.out.bam
            ch_bqsr_bai = SENTIEON_BQSR.out.bai
            ch_bqsr_csv = SENTIEON_BQSR.out.recal_csv
            ch_versions = ch_versions.mix(SENTIEON_BQSR.out.versions.first())
        }

        ch_versions = ch_versions.mix(SENTIEON_BWAMEM.out.versions.first())
        ch_versions = ch_versions.mix(SENTIEON_DATAMETRICS.out.versions.first())
        ch_versions = ch_versions.mix(SENTIEON_LOCUSCOLLECTOR.out.versions.first())
        ch_versions = ch_versions.mix(SENTIEON_DEDUP.out.versions.first())

    emit:
        marked_bam  = SENTIEON_DEDUP.out.bam                             // channel: [ val(meta), path(bam) ]
        marked_bai  = SENTIEON_DEDUP.out.bai                             // channel: [ val(meta), path(bai) ]
        recal_bam   = ch_bqsr_bam.ifEmpty(null)                          // channel: [ val(meta), path(bam) ]
        recal_bai   = ch_bqsr_bai.ifEmpty(null)                          // channel: [ val(meta), path(bai) ]
        recal_csv   = ch_bqsr_csv.ifEmpty(null)                          // channel: [ val(meta), path(csv) ]
        mq_metrics  = SENTIEON_DATAMETRICS.out.mq_metrics.ifEmpty(null)  // channel: [ val(meta), path(mq_metrics) ]
        qd_metrics  = SENTIEON_DATAMETRICS.out.qd_metrics.ifEmpty(null)  // channel: [ val(meta), path(qd_metrics) ]
        gc_metrics  = SENTIEON_DATAMETRICS.out.gc_metrics.ifEmpty(null)  // channel: [ val(meta), path(gc_metrics) ]
        gc_summary  = SENTIEON_DATAMETRICS.out.gc_summary.ifEmpty(null)  // channel: [ val(meta), path(gc_summary) ]
        aln_metrics = SENTIEON_DATAMETRICS.out.aln_metrics.ifEmpty(null) // channel: [ val(meta), path(aln_metrics) ]
        is_metrics  = SENTIEON_DATAMETRICS.out.is_metrics.ifEmpty(null)  // channel: [ val(meta), path(is_metrics) ]
        versions    = ch_versions                                        // channel: [ path(versions.yml) ]
}
