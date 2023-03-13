//
// A subworkflow to create genome interval files necessary for bam/vcf scatter operations.
//

include { BUILD_BED            } from '../../modules/local/create_bed_from_fai'
include { GATK4_SPLITINTERVALS } from '../../modules/nf-core/gatk4/splitintervals/main'

workflow SCATTER_GENOME {

    take:
        ch_dict               // channel: [mandatory] [ path(dict) ]
        ch_fai_meta           // channel: [mandatory] [ val(meta), path(fai) ]
        ch_fai_no_meta        // channel: [mandatory] [ path(fai) ]
        ch_fasta_no_meta      // channel: [mandatory] [ path(fasta) ]

    main:
        ch_versions = Channel.empty()

        BUILD_BED (ch_fai_meta)

        GATK4_SPLITINTERVALS(BUILD_BED.out.bed, ch_fasta_no_meta, ch_fai_no_meta, ch_dict)

        ch_versions = ch_versions.mix(BUILD_BED.out.versions)
        ch_versions = ch_versions.mix(GATK4_SPLITINTERVALS.out.versions)
        GATK4_SPLITINTERVALS.out.split_intervals.view()

    emit:
        bed             = BUILD_BED.out.bed.collect()   // channel: [ val(meta), path(bed) ]
        split_intervals = GATK4_SPLITINTERVALS.out.split_intervals.map { meta, it -> it }.flatten().collate(1) // channel: [ val(meta), [ path(interval_lists) ] ]
        versions        = ch_versions                   // channel: [ path(versions.yml) ]
}
