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

def tabix_vcf_check                         = params.tabix_options.clone()
tabix_vcf_check.publish_dir                 = "vcf_check/"

include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS } from '../../modules/nf-core/modules/bcftools/norm/main'  addParams( options: split_multiallelics_vcf_check )
include { BCFTOOLS_NORM as REMOVE_DUPLICATES } from '../../modules/nf-core/modules/bcftools/norm/main'  addParams( options: rm_duplicates_vcf_check )
include { CHECK_INPUT_VCF } from '../../modules/local/check_input_vcf' addParams( options: params.vcf_options )
include { TABIX_TABIX as TABIX } from '../../modules/nf-core/modules/tabix/tabix/main'  addParams( options: tabix_vcf_check )

workflow CHECK_VCF {
    take:
    vcfs   // array: [ vcf files ]
    fasta  // path(fasta)

    main:
    CHECK_INPUT_VCF( vcfs )
        .filter { it.size()>0 }
        .splitCsv()
        .map { [ [ 'id':it[0] ], it[1] ] }
        .set{ ch_unprocessed_vcfs }

    ch_unprocessed_vcfs.view()

    SPLIT_MULTIALLELICS (ch_unprocessed_vcfs, fasta)
    REMOVE_DUPLICATES (SPLIT_MULTIALLELICS.out.vcf, fasta)
    TABIX (REMOVE_DUPLICATES.out.vcf)

    emit:
        vcf  =  REMOVE_DUPLICATES.out.vcf  // path: genome.fasta
        idx  =  TABIX.out.tbi
}
