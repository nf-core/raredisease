//
// call Single-nucleotide Varinats
//

include { CALL_SNV_DEEPVARIANT             } from './variant_calling/call_snv_deepvariant'
include { CALL_SNV_SENTIEON                } from './variant_calling/call_snv_sentieon'
include { CALL_SNV_MT                      } from './variant_calling/call_snv_MT'
include { CALL_SNV_MT as CALL_SNV_MT_SHIFT } from './variant_calling/call_snv_MT'
include { POSTPROCESS_MT_CALLS             } from './variant_calling/postprocess_MT_calls'
include { GATK4_SELECTVARIANTS             } from '../../modules/nf-core/gatk4/selectvariants/main'

workflow CALL_SNV {
    take:
        ch_genome_bam_bai     // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_mt_bam_bai         // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_mtshift_bam_bai    // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_chrsizes    // channel: [mandatory] [ path(sizes) ]
        ch_genome_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai         // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_dictionary  // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mt_intervals       // channel: [optional] [ path(interval_list) ]
        ch_mtshift_fasta      // channel: [optional] [ val(meta), path(fasta) ]
        ch_mtshift_fai        // channel: [optional] [ val(meta), path(fai) ]
        ch_mtshift_dictionary // channel: [optional] [ val(meta), path(dict) ]
        ch_mtshift_intervals  // channel: [optional] [ path(interval_list) ]
        ch_mtshift_backchain  // channel: [mandatory] [ val(meta), path(back_chain) ]
        ch_dbsnp              // channel: [optional] [ val(meta), path(vcf) ]
        ch_dbsnp_tbi          // channel: [optional] [ val(meta), path(tbi) ]
        ch_call_interval      // channel: [mandatory] [ path(intervals) ]
        ch_ml_model           // channel: [mandatory] [ path(model) ]
        ch_case_info          // channel: [mandatory] [ val(case_info) ]
        ch_foundin_header     // channel: [mandatory] [ path(header) ]
        ch_pcr_indel_model    // channel: [optional] [ val(sentieon_dnascope_pcr_indel_model) ]

    main:
        ch_versions     = Channel.empty()
        ch_deepvar_vcf  = Channel.empty()
        ch_deepvar_tbi  = Channel.empty()
        ch_sentieon_vcf = Channel.empty()
        ch_sentieon_tbi = Channel.empty()

        if (params.variant_caller.equals("deepvariant")) {
            CALL_SNV_DEEPVARIANT (      // triggered only when params.variant_caller is set as deepvariant
                ch_genome_bam_bai,
                ch_genome_fasta,
                ch_genome_fai,
                ch_case_info,
                ch_foundin_header,
                ch_genome_chrsizes
            )
            ch_deepvar_vcf = CALL_SNV_DEEPVARIANT.out.vcf
            ch_deepvar_tbi = CALL_SNV_DEEPVARIANT.out.tabix
            ch_versions    = ch_versions.mix(CALL_SNV_DEEPVARIANT.out.versions)
        } else if (params.variant_caller.equals("sentieon")) {
            CALL_SNV_SENTIEON(         // triggered only when params.variant_caller is set as sentieon
                ch_genome_bam_bai,
                ch_genome_fasta,
                ch_genome_fai,
                ch_dbsnp,
                ch_dbsnp_tbi,
                ch_call_interval,
                ch_ml_model,
                ch_case_info,
                ch_pcr_indel_model,
                ch_foundin_header,
                ch_genome_chrsizes
            )
            ch_sentieon_vcf = CALL_SNV_SENTIEON.out.vcf
            ch_sentieon_tbi = CALL_SNV_SENTIEON.out.tabix
            ch_versions    = ch_versions.mix(CALL_SNV_SENTIEON.out.versions)
        }

        ch_vcf       = Channel.empty().mix(ch_deepvar_vcf, ch_sentieon_vcf)
        ch_tabix     = Channel.empty().mix(ch_deepvar_tbi, ch_sentieon_tbi)

        ch_vcf
            .join(ch_tabix, failOnMismatch:true, failOnDuplicate:true)
            .map { meta, vcf, tbi -> return [meta, vcf, tbi, []]}
            .set {ch_selvar_in}
        GATK4_SELECTVARIANTS(ch_selvar_in) // remove mitochondrial variants

        ch_genome_vcf       = GATK4_SELECTVARIANTS.out.vcf
        ch_genome_tabix     = GATK4_SELECTVARIANTS.out.tbi
        ch_genome_vcf_tabix = ch_genome_vcf.join(ch_genome_tabix, failOnMismatch:true, failOnDuplicate:true)

        CALL_SNV_MT(
            ch_mt_bam_bai,
            ch_genome_fasta,
            ch_genome_fai,
            ch_genome_dictionary,
            ch_mt_intervals
        )

        CALL_SNV_MT_SHIFT(
            ch_mtshift_bam_bai,
            ch_mtshift_fasta,
            ch_mtshift_fai,
            ch_mtshift_dictionary,
            ch_mtshift_intervals
        )

        POSTPROCESS_MT_CALLS(
            CALL_SNV_MT.out.vcf,
            CALL_SNV_MT_SHIFT.out.vcf,
            ch_genome_fasta,
            ch_genome_dictionary,
            ch_genome_fai,
            ch_mtshift_backchain,
            ch_case_info,
            ch_foundin_header,
            ch_genome_chrsizes
        )

        ch_versions = ch_versions.mix(CALL_SNV_MT.out.versions)
        ch_versions = ch_versions.mix(CALL_SNV_MT_SHIFT.out.versions)
        ch_versions = ch_versions.mix(POSTPROCESS_MT_CALLS.out.versions)
        ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions)

    emit:
        genome_vcf       = ch_genome_vcf                // channel: [ val(meta), path(vcf) ]
        genome_tabix     = ch_genome_tabix              // channel: [ val(meta), path(tbi) ]
        genome_vcf_tabix = ch_genome_vcf_tabix          // channel: [ val(meta), path(vcf), path(tbi) ]
        mt_vcf           = POSTPROCESS_MT_CALLS.out.vcf // channel: [ val(meta), path(vcf) ]
        mt_tabix         = POSTPROCESS_MT_CALLS.out.tbi // channel: [ val(meta), path(tbi) ]
        versions         = ch_versions                  // channel: [ path(versions.yml) ]
}
