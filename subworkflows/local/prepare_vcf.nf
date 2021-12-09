//
// Prepare reference vcf files
//

include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_PV } from '../../modules/nf-core/modules/bcftools/norm/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES_PV } from '../../modules/nf-core/modules/bcftools/norm/main'
include { CHECK_INPUT_VCF } from '../../modules/local/check_input_vcf'
include { TABIX_TABIX as TABIX_PV } from '../../modules/nf-core/modules/tabix/tabix/main'

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

    ch_vcfs_norm.unprocessed.view()
    SPLIT_MULTIALLELICS_PV (ch_vcfs_norm.unprocessed, fasta)

    REMOVE_DUPLICATES_PV (SPLIT_MULTIALLELICS_PV.out.vcf, fasta).vcf
        .set { ch_vcfs_rmdup }

    vcf_out = ch_vcfs_rmdup.mix( ch_vcfs_norm.processed )

    TABIX_PV (vcf_out)

    emit:
        vcf  =  vcf_out        // path: normalized_vcf
        idx  =  TABIX_PV.out.tbi
}
