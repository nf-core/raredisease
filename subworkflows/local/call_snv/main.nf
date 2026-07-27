//
// call Single-nucleotide Varinats
//

include { CALL_SNV_DEEPVARIANT } from '../call_snv_deepvariant'
include { CALL_SNV_SENTIEON    } from '../call_snv_sentieon'
include { GATK4_SELECTVARIANTS } from '../../../modules/nf-core/gatk4/selectvariants/main'

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
        ch_par_bed                // channel: [optional] [ val(meta), path(bed) ]
        ch_pcr_indel_model        // channel: [optional] [ val(sentieon_dnascope_pcr_indel_model) ]
        ch_target_bed             // channel: [mandatory] [ val(meta), path(bed), path(index) ]
        val_analysis_type             // string:  'wgs', 'wes', or 'mito'
        val_skip_split_multiallelics  // boolean
        val_variant_caller            // string:  'deepvariant' or 'sentieon'

    main:
        ch_deepvariant_gvcf     = channel.empty()
        ch_deepvariant_gtbi     = channel.empty()
        ch_deepvariant_report   = channel.empty()
        ch_deepvariant_tbi      = channel.empty()
        ch_deepvariant_vcf      = channel.empty()
        ch_sentieon_gvcf        = channel.empty()
        ch_sentieon_gtbi        = channel.empty()
        ch_sentieon_tbi         = channel.empty()
        ch_sentieon_vcf         = channel.empty()

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
                val_analysis_type,
                val_skip_split_multiallelics
            )
            ch_deepvariant_vcf    = CALL_SNV_DEEPVARIANT.out.vcf
            ch_deepvariant_tbi    = CALL_SNV_DEEPVARIANT.out.tabix
            ch_deepvariant_gvcf   = CALL_SNV_DEEPVARIANT.out.gvcf
            ch_deepvariant_gtbi   = CALL_SNV_DEEPVARIANT.out.gvcf_tabix
            ch_deepvariant_report = CALL_SNV_DEEPVARIANT.out.deepvariant_report
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
                ch_pcr_indel_model,
                val_skip_split_multiallelics
            )
            ch_sentieon_vcf  = CALL_SNV_SENTIEON.out.vcf
            ch_sentieon_tbi  = CALL_SNV_SENTIEON.out.tabix
            ch_sentieon_gvcf = CALL_SNV_SENTIEON.out.gvcf
            ch_sentieon_gtbi = CALL_SNV_SENTIEON.out.gvcf_tbi
        }

        ch_vcf    = channel.empty().mix(ch_deepvariant_vcf, ch_sentieon_vcf)
        ch_tabix  = channel.empty().mix(ch_deepvariant_tbi, ch_sentieon_tbi)
        ch_gvcf   = channel.empty().mix(ch_deepvariant_gvcf, ch_sentieon_gvcf)
        ch_gtabix = channel.empty().mix(ch_deepvariant_gtbi, ch_sentieon_gtbi)

        ch_vcf
            .join(ch_tabix, failOnMismatch:true, failOnDuplicate:true)
            .map { meta, vcf, tbi -> return [meta, vcf, tbi, []]}
            .set {ch_select_variants_in}
        GATK4_SELECTVARIANTS(ch_select_variants_in) // remove mitochondrial variants

        ch_genome_vcf       = GATK4_SELECTVARIANTS.out.vcf
        ch_genome_tabix     = GATK4_SELECTVARIANTS.out.tbi
        ch_genome_vcf_tabix = ch_genome_vcf.join(ch_genome_tabix, failOnMismatch:true, failOnDuplicate:true)

    emit:
        deepvariant_report  = ch_deepvariant_report   // channel: [ val(meta), path(html) ]
        genome_gtabix       = ch_gtabix               // channel: [ val(meta), path(gtbi) ]
        genome_gvcf         = ch_gvcf                 // channel: [ val(meta), path(gvcf) ]
        genome_tabix        = ch_genome_tabix         // channel: [ val(meta), path(tbi) ]
        genome_vcf          = ch_genome_vcf           // channel: [ val(meta), path(vcf) ]
        genome_vcf_tabix    = ch_genome_vcf_tabix     // channel: [ val(meta), path(vcf), path(tbi) ]
}
