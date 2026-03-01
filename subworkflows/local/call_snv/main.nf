//
// call Single-nucleotide Varinats
//

include { CALL_SNV_DEEPVARIANT             } from '../call_snv_deepvariant'
include { CALL_SNV_SENTIEON                } from '../call_snv_sentieon'
include { CALL_SNV_MT                      } from '../call_snv_MT'
include { CALL_SNV_MT as CALL_SNV_MT_SHIFT } from '../call_snv_MT'
include { POSTPROCESS_MT_CALLS             } from '../postprocess_MT_calls'
include { GATK4_SELECTVARIANTS             } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { BCFTOOLS_CONCAT                  } from '../../../modules/nf-core/bcftools/concat'

workflow CALL_SNV {
    take:
        ch_call_interval          // channel: [mandatory] [ path(intervals) ]
        ch_case_info              // channel: [mandatory] [ val(case_info) ]
        ch_dbsnp                  // channel: [optional] [ val(meta), path(vcf) ]
        ch_dbsnp_tbi              // channel: [optional] [ val(meta), path(tbi) ]
        ch_foundin_header         // channel: [mandatory] [ path(header) ]
        ch_genome_bam_bai         // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_chrsizes        // channel: [mandatory] [ path(sizes) ]
        ch_genome_fasta           // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai             // channel: [mandatory] [ val(meta), path(fai) ]
        ch_ml_model               // channel: [mandatory] [ path(model) ]
        ch_mt_bam_bai             // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_mt_dictionary          // channel: [optional] [ val(meta), path(dict) ]
        ch_mt_fai                 // channel: [optional] [ val(meta), path(fai) ]
        ch_mt_fasta               // channel: [optional] [ val(meta), path(fasta) ]
        ch_mt_intervals           // channel: [optional] [ path(interval_list) ]
        ch_mtshift_bam_bai        // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_mtshift_dictionary     // channel: [optional] [ val(meta), path(dict) ]
        ch_mtshift_fai            // channel: [optional] [ val(meta), path(fai) ]
        ch_mtshift_fasta          // channel: [optional] [ val(meta), path(fasta) ]
        ch_mtshift_intervals      // channel: [optional] [ path(interval_list) ]
        ch_mtshift_backchain      // channel: [mandatory] [ val(meta), path(back_chain) ]
        ch_par_bed                // channel: [optional] [ val(meta), path(bed) ]
        ch_pcr_indel_model        // channel: [optional] [ val(sentieon_dnascope_pcr_indel_model) ]
        ch_target_bed             // channel: [mandatory] [ val(meta), path(bed), path(index) ]
        val_analysis_type         // string:  'wgs', 'wes', or 'mito'
        val_concatenate_snv_calls // boolean
        val_run_mt_for_wes        // boolean
        val_variant_caller        // string:  'deepvariant' or 'sentieon'

    main:
        ch_versions      = channel.empty()
        ch_deepvar_vcf   = channel.empty()
        ch_deepvar_tbi   = channel.empty()
        ch_deepvar_gvcf  = channel.empty()
        ch_deepvar_gtbi  = channel.empty()
        ch_mt_vcf        = channel.empty()
        ch_mt_tabix      = channel.empty()
        ch_mt_vcf_tabix  = channel.empty()
        ch_mt_txt        = channel.empty()
        ch_sentieon_vcf  = channel.empty()
        ch_sentieon_tbi  = channel.empty()
        ch_sentieon_gvcf = channel.empty()
        ch_sentieon_gtbi = channel.empty()

        if (val_variant_caller.equals("deepvariant") && !val_analysis_type.equals("mito")) {
            CALL_SNV_DEEPVARIANT (
                ch_genome_bam_bai,
                ch_case_info,
                ch_foundin_header,
                ch_genome_chrsizes,
                ch_genome_fai,
                ch_genome_fasta,
                ch_par_bed,
                ch_target_bed,
                val_analysis_type
            )
            ch_deepvar_vcf  = CALL_SNV_DEEPVARIANT.out.vcf
            ch_deepvar_tbi  = CALL_SNV_DEEPVARIANT.out.tabix
            ch_deepvar_gvcf = CALL_SNV_DEEPVARIANT.out.gvcf
            ch_deepvar_gtbi = CALL_SNV_DEEPVARIANT.out.gvcf_tabix
        } else if (val_variant_caller.equals("sentieon")) {
            CALL_SNV_SENTIEON(
                ch_genome_bam_bai,
                ch_call_interval,
                ch_case_info,
                ch_dbsnp,
                ch_dbsnp_tbi,
                ch_foundin_header,
                ch_genome_chrsizes,
                ch_genome_fai,
                ch_genome_fasta,
                ch_ml_model,
                ch_pcr_indel_model
            )
            ch_sentieon_vcf  = CALL_SNV_SENTIEON.out.vcf
            ch_sentieon_tbi  = CALL_SNV_SENTIEON.out.tabix
            ch_sentieon_gvcf = CALL_SNV_SENTIEON.out.gvcf
            ch_sentieon_gtbi = CALL_SNV_SENTIEON.out.gvcf_tbi
        }

        ch_vcf    = channel.empty().mix(ch_deepvar_vcf, ch_sentieon_vcf)
        ch_tabix  = channel.empty().mix(ch_deepvar_tbi, ch_sentieon_tbi)
        ch_gvcf   = channel.empty().mix(ch_deepvar_gvcf, ch_sentieon_gvcf)
        ch_gtabix = channel.empty().mix(ch_deepvar_gtbi, ch_sentieon_gtbi)

        ch_vcf
            .join(ch_tabix, failOnMismatch:true, failOnDuplicate:true)
            .map { meta, vcf, tbi -> return [meta, vcf, tbi, []]}
            .set {ch_selvar_in}
        GATK4_SELECTVARIANTS(ch_selvar_in) // remove mitochondrial variants

        ch_genome_vcf       = GATK4_SELECTVARIANTS.out.vcf
        ch_genome_tabix     = GATK4_SELECTVARIANTS.out.tbi
        ch_genome_vcf_tabix = ch_genome_vcf.join(ch_genome_tabix, failOnMismatch:true, failOnDuplicate:true)

        if (val_analysis_type.matches("wgs|mito") || val_run_mt_for_wes) {
            CALL_SNV_MT(
                ch_mt_bam_bai,
                ch_mt_dictionary,
                ch_mt_fai,
                ch_mt_fasta,
                ch_mt_intervals
            )

            CALL_SNV_MT_SHIFT(
                ch_mtshift_bam_bai,
                ch_mtshift_dictionary,
                ch_mtshift_fai,
                ch_mtshift_fasta,
                ch_mtshift_intervals
            )

            POSTPROCESS_MT_CALLS(
                ch_case_info,
                ch_foundin_header,
                ch_genome_chrsizes,
                ch_mt_dictionary,
                ch_mt_fai,
                ch_mt_fasta,
                CALL_SNV_MT.out.vcf,
                ch_mtshift_backchain,
                CALL_SNV_MT_SHIFT.out.vcf
            )
            ch_mt_vcf       = POSTPROCESS_MT_CALLS.out.vcf
            ch_mt_tabix     = POSTPROCESS_MT_CALLS.out.tbi
            ch_mt_vcf_tabix = ch_mt_vcf.join(ch_mt_tabix, failOnMismatch:true, failOnDuplicate:true)
            ch_mt_txt       = CALL_SNV_MT.out.txt
            ch_versions     = ch_versions.mix(CALL_SNV_MT.out.versions)
            ch_versions     = ch_versions.mix(CALL_SNV_MT_SHIFT.out.versions)
        }

        if (val_concatenate_snv_calls) {
            ch_concat_vcf_in = ch_genome_vcf_tabix.concat(ch_mt_vcf_tabix).groupTuple()
            BCFTOOLS_CONCAT (
                ch_concat_vcf_in
            )
        }

    emit:
        genome_gtabix    = ch_gtabix           // channel: [ val(meta), path(gtbi) ]
        genome_gvcf      = ch_gvcf             // channel: [ val(meta), path(gvcf) ]
        genome_tabix     = ch_genome_tabix     // channel: [ val(meta), path(tbi) ]
        genome_vcf       = ch_genome_vcf       // channel: [ val(meta), path(vcf) ]
        genome_vcf_tabix = ch_genome_vcf_tabix // channel: [ val(meta), path(vcf), path(tbi) ]
        mt_tabix         = ch_mt_tabix         // channel: [ val(meta), path(tbi) ]
        mt_vcf           = ch_mt_vcf           // channel: [ val(meta), path(vcf) ]
        mt_txt           = ch_mt_txt           // channel: [ val(meta), path(txt) ]
        versions         = ch_versions         // channel: [ path(versions.yml) ]
}
