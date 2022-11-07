include { BUILD_BED            } from '../../modules/local/create_bed_from_fai'
include { GATK4_SPLITINTERVALS } from '../../modules/nf-core/gatk4/splitintervals/main'
workflow SCATTER_GENOME {

    take:
        dict
        fai_meta           // channel: [ val(meta), path(vcf) ]
        fai_no_meta
        fasta_no_meta

    main:
        ch_versions = Channel.empty()

        BUILD_BED (fai_meta)
        ch_versions = ch_versions.mix(BUILD_BED.out.versions)

        GATK4_SPLITINTERVALS(BUILD_BED.out.bed, fasta_no_meta, fai_no_meta, dict)
        ch_versions = ch_versions.mix(GATK4_SPLITINTERVALS.out.versions)

    emit:
        bed             = BUILD_BED.out.bed.collect()
        split_intervals = GATK4_SPLITINTERVALS.out.split_intervals.map { meta, it -> it }.flatten().collate(1)
        versions        = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
