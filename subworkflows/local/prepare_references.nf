//
// Prepare reference files
//

include { CHECK_BED      } from './prepare_bed'
include { CHECK_VCF      } from './prepare_vcf'
include { PREPARE_GENOME } from './prepare_genome'


workflow PREPARE_REFERENCES {
    take:

    main:
        //
        // Prepare genome
        //
        ch_versions = Channel.empty()
        PREPARE_GENOME (
            params.fasta,
            params.variant_catalog
        )
        .set { ch_genome }
        ch_versions = ch_versions.mix(ch_genome.versions)

        //
        // Gnomad vcf
        //
        ch_gnomad_vcf = Channel.empty()
        ch_gnomad_idx = Channel.empty()
        if (params.gnomad) {
            CHECK_VCF(
                params.gnomad,
                ch_genome.fasta
            )
            ch_gnomad_vcf = CHECK_VCF.out.vcf
            ch_gnomad_idx = CHECK_VCF.out.idx
        }

        //
        // Target bed
        //
        ch_target_bed       = Channel.empty()
        ch_target_intervals = Channel.empty()
        ch_bait_intervals   = Channel.empty()
        if (params.target_bed) {
            CHECK_BED(
                params.target_bed,
                ch_genome.sequence_dict
            )
            ch_target_bed       = CHECK_BED.out.bed
            ch_target_intervals = CHECK_BED.out.target_intervals
            ch_bait_intervals   = CHECK_BED.out.bait_intervals
        }

    emit:
        bwamem2_index     = ch_genome.bwamem2_index
        chrom_sizes       = ch_genome.chrom_sizes
        genome_fasta      = ch_genome.fasta
        genome_fai        = ch_genome.fai
        sequence_dict     = ch_genome.sequence_dict
        variant_catalog   = ch_genome.variant_catalog
        vcfanno_resources = ch_genome.vcfanno_resources
        gnomad_vcf        = ch_gnomad_vcf
        gnomad_idx        = ch_gnomad_idx
        target_bed        = ch_target_bed
        target_intervals  = ch_target_intervals
        bait_intervals    = ch_bait_intervals
        versions          = ch_versions
}

