//
// A preprocessing workflow for Gens
//

include { GATK4_COLLECTREADCOUNTS as COLLECTREADCOUNTS } from '../../modules/local/gatk4/collectreadcounts/main'
include { GATK4_DENOISEREADCOUNTS as DENOISEREADCOUNTS } from '../../modules/local/gatk4/denoisereadcounts/main'
include { GATK4_HAPLOTYPECALLER as HAPLOTYPECALLER } from '../../modules/nf-core/modules/gatk4/haplotypecaller/main'

workflow GENS {
    take:
        bam             // channel: [ val(meta), path(bam) ]
        bai             // channel: [ val(meta), path(bai) ]
        fasta           // path(fasta)
        fai             // path(fai)
        interval_list   // path(interval_list)
        pon             // path(pon)
        gnomad_pos      // path(gnomad_pos)
        case_info       // channel: [ val(case_info) ]
        seq_dict        // path: seq_dict

    main:
        HAPLOTYPECALLER ( bam.join(bai, by: [0]), fasta, fai, seq_dict, [], [] )

        COLLECTREADCOUNTS ( bam, fasta, interval_list )

        DENOISEREADCOUNTS ( COLLECTREADCOUNTS.out.read_counts, fasta, pon )

        GENS ( DENOISEREADCOUNTS.out.standardized_read_counts, HAPLOTYPECALLER.out.vcf, gnomad_pos )

    emit:
        gens_cov_bed_gz = GENS.out.cov
        gens_baf_bed_gz = GENS.out.baf
        versions        = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
