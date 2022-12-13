//
// call Single-nucleotide Varinats
//

include { CALL_SNV_DEEPVARIANT } from './variant_calling/call_snv_deepvariant'
include { CALL_SNV_SENTIEON    } from './variant_calling/call_snv_sentieon'
include { GATK4_SELECTVARIANTS } from '../../modules/nf-core/gatk4/selectvariants/main'


workflow CALL_SNV {
    take:
	    variant_caller       // string:  params.variant_caller
	    input                // channel: [ val(meta), path(bam), path(bai) ]
	    fasta                // channel: [genome.fasta]
	    fai                  // channel: [genome.fai]
	    known_dbsnp          // channel: [ /path/to/known_dbsnp ]
	    known_dbsnp_tbi      // channel: [ /path/to/known_dbsnp_tbi ]
        call_interval        // channel: [ /path/to/call_intervals ]
	    ml_model             // channel: [ /path/to/ml_model ]
	    case_info            // channel: [ case_id ]

    main:
        ch_versions   = Channel.empty()

        if (variant_caller == "deepvariant") {
            CALL_SNV_DEEPVARIANT ( input, fasta, fai, case_info )
            ch_vcf      = CALL_SNV_DEEPVARIANT.out.vcf
            ch_tabix    = CALL_SNV_DEEPVARIANT.out.tabix
            ch_versions = ch_versions.mix(CALL_SNV_DEEPVARIANT.out.versions)
        } else if ( variant_caller == "sentieon" ) {
            CALL_SNV_SENTIEON( input, fasta, fai, known_dbsnp, known_dbsnp_tbi, call_interval, ml_model, case_info )
            ch_vcf      = CALL_SNV_SENTIEON.out.vcf
            ch_tabix    = CALL_SNV_SENTIEON.out.tabix
            ch_versions = ch_versions.mix(CALL_SNV_SENTIEON.out.versions)
        } else {
           exit 1, 'Please provide a valid variant caller!'
        }

        ch_vcf.join(ch_tabix)
            .map { meta, vcf, tbi -> return [meta, vcf, tbi, []]}
            .set { ch_selvar_in }
        GATK4_SELECTVARIANTS(ch_selvar_in)

    emit:
        vcf      = GATK4_SELECTVARIANTS.out.vcf
        tabix    = GATK4_SELECTVARIANTS.out.tbi
        versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
