//
// Prepare bam files for MT allignment
//

include { BWAMEM2_MEM as BWAMEM2_MEM_MT                         } from '../../modules/nf-core/modules/bwamem2/mem/main'
include { GATK4_MERGEBAMALIGNMENT as GATK4_MERGEBAMALIGNMENT_MT } from '../../modules/nf-core/modules/gatk4/mergebamalignment/main'
include { PICARD_MARKDUPLICATES as PICARD_MARKDUPLICATES_MT     } from '../../modules/nf-core/modules/picard/markduplicates/main'
include { HAPLOCHECK as HAPLOCHECK_MT                           } from '../../modules/nf-core/modules/haplocheck/main'
include { GATK4_MUTECT2 as GATK4_MUTECT2_MT                     } from '../../modules/nf-core/modules/gatk4/mutect2/main'

workflow ALIGN_MT {
    take:
        fastq  // TO DO d: and file: bam index: bam.bai
        fasta
        fai
        dict

    main:
        ch_versions = Channel.empty()

        // Outputs bam files
        BWAMEM2_MEM_MT ( fastq , fasta, true)
        ch_versions = ch_versions.mix(BWAMEM2_MEM_MT.out.versions.first())

        // Merges bam files
        GATK4_MERGEBAMALIGNMENT_MT ( BWAMEM2_MEM_MT.out.bam, fasta, dict )
        ch_versions = ch_versions.mix(GATK4_MERGEBAMALIGNMENT_MT.out.versions.first())

        // Marks duplicates
        PICARD_MARKDUPLICATES_MT ( GATK4_MERGEBAMALIGNMENT_MT.out.bam )
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES_MT.out.versions.first())
        
        // Calls variants with Mutect2
        GATK4_MUTECT2_MT ( PICARD_MARKDUPLICATES_MT.out.bam, fasta, fai, dict, [], [], [], [] )
        ch_versions = ch_versions.mix(GATK4_MUTECT2_MT.out.versions.first())

        // Haplocheck
        HAPLOCHECK_MT ( GATK4_MUTECT2_MT.out.vcf )
        ch_versions = ch_versions.mix(HAPLOCHECK_MT.out.versions.first())


    emit:
        vcf      = GATK4_MUTECT2_MT.out.vcf
        tbi      = GATK4_MUTECT2_MT.out.tbi
        txt      = HAPLOCHECK_MT.out.txt
        html     = HAPLOCHECK_MT.out.html
        versions = ch_versions // channel: [ versions.yml ]
}
