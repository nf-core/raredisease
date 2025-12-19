//
// Prepare bam files for MT allignment
//

include { GATK4_PRINTREADS as GATK4_PRINTREADS_MT } from '../../../modules/nf-core/gatk4/printreads/main'
include { GATK4_REVERTSAM as GATK4_REVERTSAM_MT   } from '../../../modules/nf-core/gatk4/revertsam/main'
include { GATK4_SAMTOFASTQ as GATK4_SAMTOFASTQ_MT } from '../../../modules/nf-core/gatk4/samtofastq/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MT     } from '../../../modules/nf-core/samtools/index/main'

workflow CONVERT_MT_BAM_TO_FASTQ {
    take:
        ch_bam_bai      // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_dict  // channel: [mandatory] [ val(meta), path(dict) ]
        ch_genome_fai   // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_fasta // channel: [mandatory] [ val(meta), path(fasta) ]

    main:
        ch_versions = channel.empty()

        // Outputs bam containing only MT
        GATK4_PRINTREADS_MT ( ch_bam_bai, ch_genome_fasta, ch_genome_fai, ch_genome_dict )
        SAMTOOLS_INDEX_MT (GATK4_PRINTREADS_MT.out.bam)

        ch_bam_bai = GATK4_PRINTREADS_MT.out.bam.join(SAMTOOLS_INDEX_MT.out.bai)
        // Removes alignment information
        GATK4_REVERTSAM_MT ( GATK4_PRINTREADS_MT.out.bam )

        // Outputs fastq files
        GATK4_SAMTOFASTQ_MT ( GATK4_REVERTSAM_MT.out.bam )

        ch_versions = ch_versions.mix(GATK4_PRINTREADS_MT.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_REVERTSAM_MT.out.versions.first())
        ch_versions = ch_versions.mix(GATK4_SAMTOFASTQ_MT.out.versions.first())

    emit:
        bam_bai  = ch_bam_bai                    // channel: [ val(meta), path(bam), path(bai) ]
        fastq    = GATK4_SAMTOFASTQ_MT.out.fastq // channel: [ val(meta), [ path(fastq) ] ]
        ubam     = GATK4_REVERTSAM_MT.out.bam    // channel: [ val(meta), path(bam) ]
        versions = ch_versions                   // channel: [ path(versions.yml) ]
}
