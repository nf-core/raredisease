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
        ch_versions = Channel.empty()

        EKLIPSE(ch_bam_bai,[])

        MT_DELETION(ch_bam_bai, ch_fasta)

        ch_versions = ch_versions.mix(EKLIPSE.out.versions.first())
        ch_versions = ch_versions.mix(MT_DELETION.out.versions.first())

    emit:
        eklipse_del    = EKLIPSE.out.deletions         // channel: [ val(meta), path(csv) ]
        eklipse_genes  = EKLIPSE.out.genes             // channel: [ val(meta), path(csv) ]
        eklipse_circos = EKLIPSE.out.circos            // channel: [ val(meta), path(png) ]
        mt_del_result  = MT_DELETION.out.mt_del_result // channel: [ val(meta), path(txt) ]
        versions       = ch_versions                   // channel: [ path(versions.yml) ]
}
