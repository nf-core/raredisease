//
// A preprocessing workflow for Gens
//

include { GATK4_COLLECTREADCOUNTS as COLLECTREADCOUNTS } from '../../modules/local/gatk4/collectreadcounts/main'
include { GATK4_DENOISEREADCOUNTS as DENOISEREADCOUNTS } from '../../modules/local/gatk4/denoisereadcounts/main'
include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER } from '../../modules/nf-core/modules/gatk4/haplotypecaller/main'
include { GENS as GENS_GENERATE } from '../../modules/local/gens/main'

workflow GENS {
    take:
        bam             // channel: [ val(meta), path(bam), path(bai) ]
        fasta           // path(fasta)
        fai             // path(fai)
        interval_list   // path(interval_list)
        pon             // path(pon)
        gnomad_pos      // path(gnomad_pos)
        case_info       // channel: [ val(case_info) ]
        seq_dict        // path: seq_dict

    main:
        ch_versions = Channel.empty()
        bam.map { meta, bam, bai ->
                        return [meta, bam, bai, []]
            }
            .set { ch_bam }

        HAPLOTYPECALLER ( ch_bam, fasta, fai, seq_dict, [], [] )
        ch_versions = ch_versions.mix(HAPLOTYPECALLER.out.versions)

        COLLECTREADCOUNTS ( ch_bam, fasta, fai, seq_dict, interval_list )
        ch_versions = ch_versions.mix(COLLECTREADCOUNTS.out.versions)

        DENOISEREADCOUNTS ( COLLECTREADCOUNTS.out.read_counts, pon )
        ch_versions = ch_versions.mix(DENOISEREADCOUNTS.out.versions)

        GENS_GENERATE ( DENOISEREADCOUNTS.out.standardized_read_counts, HAPLOTYPECALLER.out.vcf, gnomad_pos )
        ch_versions = ch_versions.mix(GENS_GENERATE.out.versions)

    emit:
        gens_cov_bed_gz = GENS_GENERATE.out.cov
        gens_baf_bed_gz = GENS_GENERATE.out.baf
        versions        = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
