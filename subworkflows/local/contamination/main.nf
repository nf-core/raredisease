//
// Subworkflow: Check sample contamination using VerifyBamID2 and/or GATK
//

include { GATK4_CALCULATECONTAMINATION } from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_GETPILEUPSUMMARIES     } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { PARSE_CONTAMINATION          } from '../../../modules/local/parse_contamination/main'
include { VERIFYBAMID_VERIFYBAMID2     } from '../../../modules/nf-core/verifybamid/verifybamid2/main'

workflow CONTAMINATION {

    take:
        ch_bam_bai             // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_fasta               // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_fai                 // channel: [mandatory] [ val(meta), path(fai) ]
        ch_dict                // channel: [optional]  [ val(meta), path(dict) ]
        ch_svd_bed             // channel: [optional]  [ path(bed) ]
        ch_svd_mu              // channel: [optional]  [ path(mu) ]
        ch_svd_ud              // channel: [optional]  [ path(ud) ]
        ch_contamination_sites // channel: [optional]  [ path(vcf), path(tbi) ]
        ch_intervals           // channel: [optional]  [ path(bed) ]
        skip_gatkcontamination // boolean: [mandatory] skip GATK contamination
        skip_verifybamid       // boolean: [mandatory] skip VerifyBamID2

    main:
        ch_gatk_contamination_mqc    = channel.empty()
        ch_gatk_contamination_pileup = channel.empty()
        ch_gatk_contamination_table  = channel.empty()
        ch_verifybamid_ancestry      = channel.empty()
        ch_verifybamid_bed           = channel.empty()
        ch_verifybamid_log           = channel.empty()
        ch_verifybamid_mu            = channel.empty()
        ch_verifybamid_self_sm       = channel.empty()
        ch_verifybamid_ud            = channel.empty()

        if (!skip_verifybamid) {
            ch_svd_in = ch_svd_ud.combine(ch_svd_mu).combine(ch_svd_bed).collect()
            VERIFYBAMID_VERIFYBAMID2(ch_bam_bai, ch_svd_in, [], ch_fasta.map { _meta, fasta -> fasta })
            ch_verifybamid_ancestry = VERIFYBAMID_VERIFYBAMID2.out.ancestry
            ch_verifybamid_bed      = VERIFYBAMID_VERIFYBAMID2.out.bed
            ch_verifybamid_log      = VERIFYBAMID_VERIFYBAMID2.out.log
            ch_verifybamid_mu       = VERIFYBAMID_VERIFYBAMID2.out.mu
            ch_verifybamid_self_sm  = VERIFYBAMID_VERIFYBAMID2.out.self_sm
            ch_verifybamid_ud       = VERIFYBAMID_VERIFYBAMID2.out.ud
        }

        if (!skip_gatkcontamination) {
            // Prepare BAM with intervals - conditionally based on analysis type
            // For WGS: intervals will be empty [], for WES: intervals will contain the BED file
            ch_bam_with_intervals = ch_bam_bai
                .combine(ch_intervals.ifEmpty([[]]))
                .map { meta, bam, bai, bed ->
                    // If bed is an empty list, pass empty list; otherwise pass bed file
                    def intervals = (bed instanceof List && bed.isEmpty()) ? [] : (bed ?: [])
                    [ meta, bam, bai, intervals ]
                }

            ch_variants_vcf = ch_contamination_sites.map { vcf, tbi -> vcf }.collect()
            ch_variants_tbi = ch_contamination_sites.map { vcf, tbi -> tbi }.collect()

            GATK4_GETPILEUPSUMMARIES(
                ch_bam_with_intervals,
                ch_fasta,
                ch_fai,
                ch_dict,
                ch_variants_vcf,
                ch_variants_tbi
            )

            GATK4_CALCULATECONTAMINATION(
                GATK4_GETPILEUPSUMMARIES.out.table.map { meta, table -> [ meta, table, [] ] }
            )

            PARSE_CONTAMINATION(GATK4_CALCULATECONTAMINATION.out.contamination)

            ch_gatk_contamination_mqc    = PARSE_CONTAMINATION.out.mqc_table
            ch_gatk_contamination_pileup = GATK4_GETPILEUPSUMMARIES.out.table
            ch_gatk_contamination_table  = GATK4_CALCULATECONTAMINATION.out.contamination
        }

    emit:
        gatk_contamination_mqc    = ch_gatk_contamination_mqc    // channel: [ val(meta), path(tsv) ]
        gatk_contamination_pileup = ch_gatk_contamination_pileup // channel: [ val(meta), path(table) ]
        gatk_contamination_table  = ch_gatk_contamination_table  // channel: [ val(meta), path(table) ]
        verifybamid_ancestry      = ch_verifybamid_ancestry      // channel: [ val(meta), path(ancestry) ]
        verifybamid_bed           = ch_verifybamid_bed           // channel: [ val(meta), path(bed) ]
        verifybamid_log           = ch_verifybamid_log           // channel: [ val(meta), path(log) ]
        verifybamid_mu            = ch_verifybamid_mu            // channel: [ val(meta), path(mu) ]
        verifybamid_self_sm       = ch_verifybamid_self_sm       // channel: [ val(meta), path(selfSM) ]
        verifybamid_ud            = ch_verifybamid_ud            // channel: [ val(meta), path(ud) ]
}
