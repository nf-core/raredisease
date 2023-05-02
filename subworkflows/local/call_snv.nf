//
// call Single-nucleotide Varinats
//

include { CALL_SNV_DEEPVARIANT } from './variant_calling/call_snv_deepvariant'
include { CALL_SNV_SENTIEON    } from './variant_calling/call_snv_sentieon'
include { GATK4_SELECTVARIANTS } from '../../modules/nf-core/gatk4/selectvariants/main'


workflow CALL_SNV {
    take:
        ch_input           // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_fasta           // channel: [mandatory] [ path(fasta) ]
        ch_fai             // channel: [mandatory] [ path(fai) ]
        ch_known_dbsnp     // channel: [optional] [ val(meta), path(vcf) ]
        ch_known_dbsnp_tbi // channel: [optional] [ val(meta), path(tbi) ]
        ch_call_interval   // channel: [mandatory] [ path(intervals) ]
        ch_ml_model        // channel: [mandatory] [ path(model) ]
        ch_case_info       // channel: [mandatory] [ val(case_info) ]

    main:
        ch_versions   = Channel.empty()
        ch_vcf        = Channel.empty()
        ch_tabix      = Channel.empty()

        CALL_SNV_DEEPVARIANT (
            ch_input,
            ch_fasta,
            ch_fai,
            ch_case_info
        )

        CALL_SNV_SENTIEON(
            ch_input,
            ch_fasta,
            ch_fai,
            ch_known_dbsnp,
            ch_known_dbsnp_tbi,
            ch_call_interval,
            ch_ml_model,
            ch_case_info
        )

        ch_vcf      = Channel.empty().mix(CALL_SNV_DEEPVARIANT.out.vcf, CALL_SNV_SENTIEON.out.vcf)
        ch_tabix    = Channel.empty().mix(CALL_SNV_DEEPVARIANT.out.tabix, CALL_SNV_SENTIEON.out.tabix)

        ch_versions = ch_versions.mix(CALL_SNV_DEEPVARIANT.out.versions)
        ch_versions = ch_versions.mix(CALL_SNV_SENTIEON.out.versions)

    emit:
        vcf      = ch_vcf      // channel: [ val(meta), path(vcf) ]
        tabix    = ch_tabix    // channel: [ val(meta), path(tbi) ]
        versions = ch_versions // channel: [ path(versions.yml) ]
}
