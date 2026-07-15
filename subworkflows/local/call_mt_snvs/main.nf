//
// call mitochondrial SNVs
//

include { CALL_SNV_MT                      } from '../call_snv_MT'
include { CALL_SNV_MT as CALL_SNV_MT_SHIFT } from '../call_snv_MT'
include { POSTPROCESS_MT_CALLS             } from '../postprocess_MT_calls'

workflow CALL_MT_SNVS {
    take:
        ch_case_info          // channel: [mandatory] [ val(case_info) ]
        ch_foundin_header     // channel: [mandatory] [ path(header) ]
        ch_genome_chrsizes    // channel: [mandatory] [ path(sizes) ]
        ch_mt_bam_bai         // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_mt_dictionary      // channel: [optional] [ val(meta), path(dict) ]
        ch_mt_fai             // channel: [optional] [ val(meta), path(fai) ]
        ch_mt_fasta           // channel: [optional] [ val(meta), path(fasta) ]
        ch_mt_intervals       // channel: [optional] [ path(interval_list) ]
        ch_mtshift_backchain  // channel: [mandatory] [ val(meta), path(back_chain) ]
        ch_mtshift_bam_bai    // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_mtshift_dictionary // channel: [optional] [ val(meta), path(dict) ]
        ch_mtshift_fai        // channel: [optional] [ val(meta), path(fai) ]
        ch_mtshift_fasta      // channel: [optional] [ val(meta), path(fasta) ]
        ch_mtshift_intervals  // channel: [optional] [ path(interval_list) ]

    main:
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

        ch_vcf     = POSTPROCESS_MT_CALLS.out.vcf
        ch_tbi     = POSTPROCESS_MT_CALLS.out.tbi
        ch_vcf_tbi = ch_vcf.join(ch_tbi, failOnMismatch:true, failOnDuplicate:true)

    emit:
        vcf     = ch_vcf      // channel: [ val(meta), path(vcf) ]
        tbi     = ch_tbi      // channel: [ val(meta), path(tbi) ]
        vcf_tbi = ch_vcf_tbi  // channel: [ val(meta), path(vcf), path(tbi) ]
}
