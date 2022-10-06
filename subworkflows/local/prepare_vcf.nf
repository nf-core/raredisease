//
// Prepare reference vcf files
//

include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_PV } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES_PV   } from '../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX_PV                 } from '../../modules/nf-core/tabix/tabix/main'
include { CHECK_INPUT_VCF                         } from '../../modules/local/check_input_vcf'

workflow CHECK_VCF {
    take:
        vcf    // file: vcf file
        fasta  // path(fasta)

    main:
        vcf_file = file(vcf)
        ch_versions = Channel.empty()

        CHECK_INPUT_VCF( vcf_file )
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

        SPLIT_MULTIALLELICS_PV (ch_vcfs_norm.unprocessed, fasta)
        ch_versions = ch_versions.mix(SPLIT_MULTIALLELICS_PV.out.versions)

        REMOVE_DUPLICATES_PV (SPLIT_MULTIALLELICS_PV.out.vcf, fasta)
            .vcf
            .set { ch_vcfs_rmdup }
        ch_versions = ch_versions.mix(REMOVE_DUPLICATES_PV.out.versions)

        vcf_out = ch_vcfs_rmdup.mix( ch_vcfs_norm.processed )

        TABIX_PV (vcf_out)

    emit:
        vcf      =  vcf_out        // path: normalized_vcf
        idx      =  TABIX_PV.out.tbi
        versions =  ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
