//
// Prepare reference vcf files
//

params.split_multiallelics_options = [:]
params.rm_duplicates_options = [:]
params.vcf_options = [:]
params.tabix_options = [:]

def split_multiallelics_vcf_check           = params.split_multiallelics_options.clone()
split_multiallelics_vcf_check.publish_files = "false"

def rm_duplicates_vcf_check                 = params.rm_duplicates_options.clone()
rm_duplicates_vcf_check.publish_dir         = "vcf_check/"
rm_duplicates_vcf_check.suffix              = "_split_rmdup"

def tabix_vcf_check                         = params.tabix_options.clone()
tabix_vcf_check.publish_dir                 = "vcf_check/"

include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS } from '../../modules/nf-core/modules/bcftools/norm/main'  addParams( options: split_multiallelics_vcf_check )
include { BCFTOOLS_NORM as REMOVE_DUPLICATES } from '../../modules/nf-core/modules/bcftools/norm/main'  addParams( options: rm_duplicates_vcf_check )
include { CHECK_INPUT_VCF } from '../../modules/local/check_input_vcf' addParams( options: params.vcf_options )

workflow CHECK_VCF {
    take:
    vcfs // array: [ vcf files ]

    main:
    ch_out = CHECK_INPUT_VCF( vcfs ).txt


  //  emit:
 //       txt                       = ch_out                  // path: genome.fasta
}
