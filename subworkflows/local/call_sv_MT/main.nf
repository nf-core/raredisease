//
// Call SV MT
//

include { MT_DELETION } from '../../../modules/local/mt_deletion_script'
include { EKLIPSE     } from '../../../modules/nf-core/eklipse/main'

workflow CALL_SV_MT {
    take:
        ch_bam_bai      // channel: [mandatory] [ val(meta), path(bam) ]
        ch_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]

    main:
        ch_versions       = Channel.empty()
        ch_eklipse_del    = Channel.empty()
        ch_eklipse_genes  = Channel.empty()
        ch_eklipse_circos = Channel.empty()

        if (!(params.skip_tools && params.skip_tools.split(',').contains('eklipse'))) {
            EKLIPSE(ch_bam_bai,[])
            ch_eklipse_del    = EKLIPSE.out.deletions
            ch_eklipse_genes  = EKLIPSE.out.genes
            ch_eklipse_circos = EKLIPSE.out.circos
            ch_versions = ch_versions.mix(EKLIPSE.out.versions.first())
        }
        MT_DELETION(ch_bam_bai, ch_fasta)

        ch_versions = ch_versions.mix(MT_DELETION.out.versions.first())

    emit:
        eklipse_del    = ch_eklipse_del                // channel: [ val(meta), path(csv) ]
        eklipse_genes  = ch_eklipse_genes              // channel: [ val(meta), path(csv) ]
        eklipse_circos = ch_eklipse_circos             // channel: [ val(meta), path(png) ]
        mt_del_result  = MT_DELETION.out.mt_del_result // channel: [ val(meta), path(txt) ]
        versions       = ch_versions                   // channel: [ path(versions.yml) ]
}
