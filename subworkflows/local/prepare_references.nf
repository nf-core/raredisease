//
// Prepare reference files
//

include { CHECK_BED                      } from './prepare_bed'
include { CHECK_VCF                      } from './prepare_vcf'
include { PREPARE_GENOME                 } from './prepare_genome'
include { TABIX_TABIX as TABIX_DBSNP     } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GNOMAD_AF } from '../../modules/nf-core/tabix/tabix/main'


workflow PREPARE_REFERENCES {
    take:
        aligner             // [mandatory] params.aligner
        bwa_index
        bwamem2_index       // [mandatory] bwamem2_index
        fasta               // [mandatory] genome.fasta
        fai                 // [mandatory] genome.fai
        gnomad
        gnomad_af
        gnomad_af_tbi
        known_dbsnp
        known_dbsnp_tbi
        target_bed
        variant_catalog     // [optional] variant_catalog.json
        vcfanno_resources   // [mandatory] vcfanno resource file

    main:
        // Prepare genome
        ch_versions = Channel.empty()
        PREPARE_GENOME (
            aligner,
            bwa_index,
            bwamem2_index,
            fasta,
            fai,
            variant_catalog,
            vcfanno_resources,
            true
        )
        .set { ch_genome }
        ch_versions = ch_versions.mix(ch_genome.versions)

        // Dbsnp vcf
        ch_dbsnp_vcf = Channel.empty()
        ch_dbsnp_tbi = Channel.empty()
        if (!known_dbsnp_tbi && known_dbsnp) {
            TABIX_DBSNP([[id:'dbsnp'], file(known_dbsnp)])
            ch_dbsnp_vcf = file(known_dbsnp)
            ch_dbsnp_tbi = TABIX_DBSNP.out.tbi.collect {it[1]}
            ch_versions  = ch_versions.mix(TABIX_DBSNP.out.versions)
        } else if (known_dbsnp_tbi && known_dbsnp) {
            ch_dbsnp_vcf = file(known_dbsnp)
            ch_dbsnp_tbi = file(known_dbsnp_tbi)
        }

        // Gnomad vcf
        ch_gnomad_vcf = Channel.empty()
        ch_gnomad_idx = Channel.empty()
        if (gnomad) {
            CHECK_VCF(
                gnomad,
                ch_genome.fasta
            )
            ch_gnomad_vcf = CHECK_VCF.out.vcf
            ch_gnomad_idx = CHECK_VCF.out.idx
            ch_versions   = ch_versions.mix(CHECK_VCF.out.versions)
        }

        // Gnomad tab
        ch_gnomad_af = Channel.empty()
        if (!gnomad_af_tbi && gnomad_af) {
            ch_gnomad_tab = [[id:'gnomad'], file(gnomad_af)]
            TABIX_GNOMAD_AF(ch_gnomad_tab)
            ch_gnomad_tbi = TABIX_GNOMAD_AF.out.tbi
            ch_gnomad_af = ch_gnomad_tab.join(ch_gnomad_tbi).collect{ it -> return [it[1], it[2]]}
            ch_versions   = ch_versions.mix(CHECK_VCF.out.versions)
        } else if (gnomad_af_tbi && gnomad_af) {
            ch_gnomad_af = [file(gnomad_af), file(gnomad_af_tbi)]
        }

        // Target bed
        ch_target_bed       = Channel.empty()
        ch_target_intervals = Channel.empty()
        ch_bait_intervals   = Channel.empty()
        if (target_bed) {
            CHECK_BED(
                target_bed,
                ch_genome.sequence_dict
            )
            ch_target_bed       = CHECK_BED.out.bed
            ch_target_intervals = CHECK_BED.out.target_intervals
            ch_bait_intervals   = CHECK_BED.out.bait_intervals
            ch_versions         = ch_versions.mix(CHECK_BED.out.versions)
        }

    emit:
        bwa_index         = ch_genome.bwa_index
        bwamem2_index     = ch_genome.bwamem2_index
        chrom_sizes       = ch_genome.chrom_sizes
        genome_fasta      = ch_genome.fasta
        genome_fai        = ch_genome.fai
        sequence_dict     = ch_genome.sequence_dict
        variant_catalog   = ch_genome.variant_catalog
        vcfanno_resources = ch_genome.vcfanno_resources
        known_dbsnp       = ch_dbsnp_vcf
        known_dbsnp_tbi   = ch_dbsnp_tbi
        gnomad_vcf        = ch_gnomad_vcf
        gnomad_idx        = ch_gnomad_idx
        gnomad_af         = ch_gnomad_af
        target_bed        = ch_target_bed
        target_intervals  = ch_target_intervals
        bait_intervals    = ch_bait_intervals
        versions          = ch_versions
}

