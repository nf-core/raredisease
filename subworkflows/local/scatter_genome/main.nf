//
// A subworkflow to create genome interval files necessary for bam/vcf scatter operations.
//

include { GATK4_SPLITINTERVALS      } from '../../../modules/nf-core/gatk4/splitintervals/main'
include { GAWK as GENOME_FAI_TO_BED } from '../../../modules/nf-core/gawk'

workflow SCATTER_GENOME {

    take:
        ch_genome_dictionary   // channel: [mandatory] [ val(meta), path(dict) ]
        ch_genome_fai          // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]

    main:

        GENOME_FAI_TO_BED (ch_genome_fai, [], false)

        GATK4_SPLITINTERVALS(GENOME_FAI_TO_BED.out.output, ch_genome_fasta, ch_genome_fai, ch_genome_dictionary)

    emit:
        gatk4_splitintervals_split_intervals = GATK4_SPLITINTERVALS.out.split_intervals // channel: [ val(meta), path(interval_list) ]
        genome_fai_to_bed_output             = GENOME_FAI_TO_BED.out.output.collect()  // channel: [ val(meta), path(bed) ]
}
