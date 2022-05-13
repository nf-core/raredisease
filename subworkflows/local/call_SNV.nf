//
// call Single-nucleotide Varinats
//

include { CALL_SNV_DEEPVARIANT } from '../nf-core/call_snv_deepvariant'
include { CALL_SNV_SENTIEON } from './calling_sentieon'
    
workflow CALL_SNV {
    take:
	     variant_caller      // string:  params.variant_caller
	    input                // channel: [ val(meta), path(bam), path(bai) ]
	    fasta                // channel: [genome.fasta]
	    fai                  // channel: [genome.fai]
	    known_dbsnp          // channel: [ /path/to/known_dbsnp ]
	    known_dbsnp_tbi      // channel: [ /path/to/known_dbsnp_tbi ]
	    ml_model             // channel: [ /path/to/ml_model ]
	    case_info            // channel: [ case_id ]
  
    main:
       ch_versions   = Channel.empty()

       if (variant_caller == "deepvariat") {
            CALL_SNV_DEEPVARIANT ( input, fasta, fai, case_info )
            ch_vcf       = CALL_SNV_DEEPVARIANT.out.vcf
            ch_vcf_index = CALL_SNV_DEEPVARIANT.out.tabix
            ch_versions  = ch_versions.mix(CALL_SNV_DEEPVARIANT.out.versions)
  
       } else if ( variant_caller == "sentieon" ) {
            CALL_SNV_SENTIEON( input, fasta, fai, known_dbsnp, known_dbsnp_tbi, ml_model )
            ch_vcf       = CALL_SNV_SENTIEON.out.vcf
            ch_vcf_index = CALL_SNV_SENTIEON.out.vcf_index
            ch_versions  = ch_versions.mix(CALL_SNV_SENTIEON.out.versions)
       } else {
            exit 1, 'Please provide a valid variant caller!'
       }
   
    emit:
     vcf            = ch_vcf          
	vcf_index      = ch_vcf_index            
     versions       = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
