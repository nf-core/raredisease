//
// A subworkflow to create genome interval files necessary for bam/vcf scatter operations.
//

include { GATK4_SPLITINTERVALS } from '../../../modules/nf-core/gatk4/splitintervals/main'
include { GAWK                 } from '../../../modules/nf-core/gawk'

workflow SCATTER_GENOME {

    take:
        ch_genome_dictionary   // channel: [mandatory] [ val(meta), path(dict) ]
        ch_genome_fai          // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
        val_save_reference     // bool

    main:

        GAWK (ch_genome_fai, [], false)

        GATK4_SPLITINTERVALS(GAWK.out.output, ch_genome_fasta, ch_genome_fai, ch_genome_dictionary)

        ch_split_intervals_publish = val_save_reference ? GATK4_SPLITINTERVALS.out.split_intervals : channel.empty()

    emit:
        bed                     = GAWK.out.output.collect()                                                         // channel: [ val(meta), path(bed) ]
        split_intervals         = GATK4_SPLITINTERVALS.out.split_intervals.map { _meta, it -> it }.flatten().collate(1) // channel: [ path(interval_list) ]
        split_intervals_publish = ch_split_intervals_publish                                                        // channel: [ val(meta), path(intervals) ]
}
