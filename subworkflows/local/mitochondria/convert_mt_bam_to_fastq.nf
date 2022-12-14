//
// Prepare bam files for MT allignment
//

include { GATK4_PRINTREADS as GATK4_PRINTREADS_MT } from '../../../modules/nf-core/gatk4/printreads/main'
include { GATK4_REVERTSAM as GATK4_REVERTSAM_MT   } from '../../../modules/nf-core/gatk4/revertsam/main'
include { GATK4_SAMTOFASTQ as GATK4_SAMTOFASTQ_MT } from '../../../modules/nf-core/gatk4/samtofastq/main'

workflow CONVERT_MT_BAM_TO_FASTQ {
    take:
        bam                    // channel: [ val(meta), file(bam), file(bai) ]
        genome_fasta_meta      // channel: [ [], genome.fasta ]
        genome_fai             // channel: [ genome.fai ]
        genome_dict            // channel: [ genome.dict ]

    main:
        ch_versions = Channel.empty()

        // Outputs bam containing only MT
        GATK4_PRINTREADS_MT ( bam, genome_fasta_meta, genome_fai, genome_dict )
        ch_versions = ch_versions.mix(GATK4_PRINTREADS_MT.out.versions.first())

        // Removes alignment information
        GATK4_REVERTSAM_MT ( GATK4_PRINTREADS_MT.out.bam )
        ch_versions = ch_versions.mix(GATK4_REVERTSAM_MT.out.versions.first())

        // Outputs fastq files
        GATK4_SAMTOFASTQ_MT ( GATK4_REVERTSAM_MT.out.bam )
        ch_versions = ch_versions.mix(GATK4_SAMTOFASTQ_MT.out.versions.first())

    emit:
        fastq    = GATK4_SAMTOFASTQ_MT.out.fastq
        bam      = GATK4_REVERTSAM_MT.out.bam
        versions = ch_versions // channel: [ versions.yml ]
}
