//
// A subworkflow to annotate structural variants.
//

include { SENTIEON_BWAMEM                            } from '../../modules/local/sentieon/bwamem'
include { SENTIEON_DRIVER as SENTIEON_DATAMETRICS    } from '../../modules/local/sentieon/driver'
include { SENTIEON_DRIVER as SENTIEON_LOCUSCOLLECTOR } from '../../modules/local/sentieon/driver'
include { SENTIEON_DRIVER as SENTIEON_DEDUP          } from '../../modules/local/sentieon/driver'
include { SENTIEON_DRIVER as SENTIEON_BQSR           } from '../../modules/local/sentieon/driver'

workflow ALIGN_SENTIEON {
    take:
        reads_input  // channel: [ val(meta), reads_input ]
        fasta        // path: genome.fasta
        fai          // path: genome.fai
        index        // channel: [ /path/to/bwamem2/index/ ]
        known_dbsnp  // path: params.known_dbsnp
        known_indels // path: params.known_indels
        known_mills  // path: params.known_mills

    main:
        ch_versions = Channel.empty()

        SENTIEON_BWAMEM ( reads_input, fasta, fai, index )
        ch_versions = ch_versions.mix(SENTIEON_BWAMEM.out.versions)

        SENTIEON_BWAMEM.out
            .bam
            .join(SENTIEON_BWAMEM.out.bai)
            .map { it -> it + [ [], [], [], [] ] }
            .set { ch_bam_bai }

        SENTIEON_DATAMETRICS (ch_bam_bai, fasta, fai, [], [], [] )
        ch_versions = ch_versions.mix(SENTIEON_DATAMETRICS.out.versions)

        SENTIEON_LOCUSCOLLECTOR (ch_bam_bai, fasta, fai, [], [], [] )

        ch_bam_bai
            .map { meta, bam, bai, score, score_idx, recal_pre, recal_post -> [ meta, bam, bai ] }
            .join(SENTIEON_LOCUSCOLLECTOR.out.score)
            .join(SENTIEON_LOCUSCOLLECTOR.out.score_idx)
            .map { it -> it + [ [], [] ] }
            .set { ch_bam_bai_score }

        SENTIEON_DEDUP ( ch_bam_bai_score, fasta, fai, [], [], [] )

        SENTIEON_DEDUP.out.bam
            .join(SENTIEON_DEDUP.out.bai)
            .map { it -> it + [ [], [], [], [] ] }
            .set { ch_dedup_bam_bai }

        SENTIEON_BQSR ( ch_dedup_bam_bai, fasta, fai, known_dbsnp, known_mills, known_indels )

    emit:
        bam                    = SENTIEON_DEDUP.out.bam
        bai                    = SENTIEON_DEDUP.out.bai
        marked_bam_bai         = ch_dedup_bam_bai
        recal_pre              = SENTIEON_BQSR.out.recal_pre
        mq_metrics             = SENTIEON_DATAMETRICS.out.mq_metrics.ifEmpty(null)
        qd_metrics             = SENTIEON_DATAMETRICS.out.qd_metrics.ifEmpty(null)
        gc_metrics             = SENTIEON_DATAMETRICS.out.gc_metrics.ifEmpty(null)
        gc_summary             = SENTIEON_DATAMETRICS.out.gc_summary.ifEmpty(null)
        aln_metrics            = SENTIEON_DATAMETRICS.out.aln_metrics.ifEmpty(null)
        is_metrics             = SENTIEON_DATAMETRICS.out.is_metrics.ifEmpty(null)
        versions               = ch_versions.ifEmpty(null)                           // channel: [ versions.yml ]
}
