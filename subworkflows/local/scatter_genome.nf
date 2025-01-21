//
// A subworkflow to create genome interval files necessary for bam/vcf scatter operations.
//

include { BUILD_BED            } from '../../modules/local/create_bed_from_fai'
include { GATK4_SPLITINTERVALS } from '../../modules/nf-core/gatk4/splitintervals/main'

workflow SCATTER_GENOME {

    take:
        ch_genome_dictionary   // channel: [mandatory] [ val(meta), path(dict) ]
        ch_genome_fai    // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_fasta  // channel: [mandatory] [ val(meta), path(fasta) ]

    main:
        ch_versions = Channel.empty()

        BUILD_BED (ch_genome_fai)

        GATK4_SPLITINTERVALS(BUILD_BED.out.bed, ch_genome_fasta, ch_genome_fai, ch_genome_dictionary)

        ch_versions = ch_versions.mix(BUILD_BED.out.versions)
        ch_versions = ch_versions.mix(GATK4_SPLITINTERVALS.out.versions)

    emit:
        bed             = BUILD_BED.out.bed.collect()   // channel: [ val(meta), path(bed) ]
        split_intervals = GATK4_SPLITINTERVALS.out.split_intervals.map { meta, it -> it }.flatten().collate(1) // channel: [ val(meta), [ path(interval_lists) ] ]
        versions        = ch_versions                   // channel: [ path(versions.yml) ]
}
