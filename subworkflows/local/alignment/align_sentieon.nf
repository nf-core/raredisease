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
        ch_fasta           // channel: [mandatory] [ path(genome.fasta) ]
        ch_fai             // channel: [mandatory] [ path(fai) ]
        ch_index           // channel: [mandatory] [ path(bwamem2_index) ]
        ch_known_dbsnp     // channel: [mandatory] [ path(known_dbsnp) ]
        ch_known_dbsnp_tbi // channel: [mandatory] [ path(known_dbsnp_tbi) ]
        val_platform       // string:  [mandatory] val(platform)

    main:
        ch_versions = Channel.empty()
        ch_bqsr_bam = Channel.empty()
        ch_bqsr_bai = Channel.empty()
        ch_bqsr_csv = Channel.empty()

        SENTIEON_BWAMEM ( ch_reads_input, ch_fasta, ch_fai, ch_index )

        SENTIEON_BWAMEM.out
            .bam
            .join(SENTIEON_BWAMEM.out.bai)
            .map{ meta, bam, bai ->
                new_meta            = meta.clone()
                new_meta.id         = new_meta.id.split('_')[0]
                new_meta.read_group = "\'@RG\\tID:" + new_meta.id + "\\tPL:" + val_platform + "\\tSM:" + new_meta.id + "\'"
                [new_meta, bam, bai]
                }
            .groupTuple(by: 0)   // group them bam paths with the same [ [samplename], [bam path, bam path, ..] ]
            .branch{
                single: it[1].size() == 1
                multiple: it[1].size() > 1
                }
            .set{ merge_bams_in }

        SENTIEON_READWRITER (merge_bams_in.multiple)
        ch_bam_bai = merge_bams_in.single.mix(SENTIEON_READWRITER.out.bam_bai)

        SENTIEON_DATAMETRICS (ch_bam_bai, ch_fasta, ch_fai )

        SENTIEON_LOCUSCOLLECTOR ( ch_bam_bai )

        ch_bam_bai
            .join(SENTIEON_LOCUSCOLLECTOR.out.score)
            .join(SENTIEON_LOCUSCOLLECTOR.out.score_idx)
            .set { ch_bam_bai_score }

        SENTIEON_DEDUP ( ch_bam_bai_score, ch_fasta, ch_fai )

        if (params.variant_caller == "sentieon") {
            SENTIEON_DEDUP.out.bam
                .join(SENTIEON_DEDUP.out.bai)
                .set { ch_dedup_bam_bai }
            SENTIEON_BQSR ( ch_dedup_bam_bai, ch_fasta, ch_fai, ch_known_dbsnp, ch_known_dbsnp_tbi )
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
        marked_bam             = SENTIEON_DEDUP.out.bam
        marked_bai             = SENTIEON_DEDUP.out.bai
        recal_bam              = ch_bqsr_bam.ifEmpty(null)
        recal_bai              = ch_bqsr_bai.ifEmpty(null)
        recal_csv              = ch_bqsr_csv.ifEmpty(null)
        mq_metrics             = SENTIEON_DATAMETRICS.out.mq_metrics.ifEmpty(null)
        qd_metrics             = SENTIEON_DATAMETRICS.out.qd_metrics.ifEmpty(null)
        gc_metrics             = SENTIEON_DATAMETRICS.out.gc_metrics.ifEmpty(null)
        gc_summary             = SENTIEON_DATAMETRICS.out.gc_summary.ifEmpty(null)
        aln_metrics            = SENTIEON_DATAMETRICS.out.aln_metrics.ifEmpty(null)
        is_metrics             = SENTIEON_DATAMETRICS.out.is_metrics.ifEmpty(null)
        versions               = ch_versions
}
