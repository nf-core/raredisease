//
// Prepare reference vcf files
//

params.split_multiallelics_options = [:]
params.rm_duplicates_options = [:]
params.vcf_options = [:]
params.tabix_options = [:]

def split_multiallelics_vcf_check           = params.split_multiallelics_options.clone()
split_multiallelics_vcf_check.publish_files = "false"
split_multiallelics_vcf_check.suffix        = "_split"

def rm_duplicates_vcf_check                 = params.rm_duplicates_options.clone()
rm_duplicates_vcf_check.publish_dir         = "vcf_check/"
rm_duplicates_vcf_check.suffix              = "_split_rmdup"

def tabix_vcf_check                         = params.tabix_options.clone()
tabix_vcf_check.publish_dir                 = "vcf_check/"

include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS } from '../../modules/nf-core/modules/bcftools/norm/main'  addParams( options: split_multiallelics_vcf_check )
include { BCFTOOLS_NORM as REMOVE_DUPLICATES } from '../../modules/nf-core/modules/bcftools/norm/main'  addParams( options: rm_duplicates_vcf_check )
include { CHECK_INPUT_VCF } from '../../modules/local/check_input_vcf' addParams( options: params.vcf_options )
include { TABIX_TABIX as TABIX } from '../../modules/nf-core/modules/tabix/tabix/main'  addParams( options: tabix_vcf_check )

workflow CHECK_VCF {
    take:
    vcf    // channel: [ vcf file ]
    fasta  // path(fasta)

    main:
    CHECK_INPUT_VCF( vcf )
        .splitCsv( header:true )
        .map { row ->
            def id        = "${row.id}"
            def filepath  = "${row.filepath}"
            def processed = "${row.processed}"
            tuple(id,filepath,processed)
        }
        .branch { id, filepath, processed ->
            processed: processed == 'yes'
                return [['id':id],filepath]
            unprocessed: processed == 'no'
                return [['id':id],filepath]
        }
        .set { ch_vcfs_norm }

    SPLIT_MULTIALLELICS (ch_vcfs_norm.unprocessed, fasta)

    REMOVE_DUPLICATES (SPLIT_MULTIALLELICS.out.vcf, fasta).vcf
        .set { ch_vcfs_rmdup }

    vcf_out = ch_vcfs_rmdup.mix( ch_vcfs_norm.processed )

    TABIX (vcf_out)

    emit:
        vcf  =  vcf_out        // path: normalized_vcf
        idx  =  TABIX.out.tbi
}
