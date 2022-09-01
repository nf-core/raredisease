//
// Prepare bam files for MT allignment
//

include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_MT       } from '../../modules/nf-core/modules/samtools/view/main'
include { GATK4_REVERTSAM as GATK4_REVERTSAM_MT   } from '../../modules/nf-core/modules/gatk4/revertsam/main'
include { GATK4_SAMTOFASTQ as GATK4_SAMTOFASTQ_MT } from '../../modules/nf-core/modules/gatk4/samtofastq/main'

workflow PREPARE_MT_ALIGNMENT {
    take:
        bam  // channel: [ val(meta), file(bam), file(bai) ]

    main:
        ch_versions = Channel.empty()

        // Outputs bam containing only MT
        SAMTOOLS_VIEW_MT ( bam, [] )
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW_MT.out.versions.first())

        // Removes alignment information
        GATK4_REVERTSAM_MT ( SAMTOOLS_VIEW_MT.out.bam )
        ch_versions = ch_versions.mix(GATK4_REVERTSAM_MT.out.versions.first())

        // Outputs fastq files
        GATK4_SAMTOFASTQ_MT ( GATK4_REVERTSAM_MT.out.bam )
        ch_versions = ch_versions.mix(GATK4_SAMTOFASTQ_MT.out.versions.first())

    emit:
        fastq    = GATK4_SAMTOFASTQ_MT.out.fastq
        bam      = GATK4_REVERTSAM_MT.out.bam
        versions = ch_versions // channel: [ versions.yml ]
}
