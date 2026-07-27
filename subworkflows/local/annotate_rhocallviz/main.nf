//
// A subworkflow to plot binned zygosity and RHO-regions.
//

include { BCFTOOLS_ROH                              } from '../../../modules/nf-core/bcftools/roh/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_RHOCALL    } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_UNCOMPRESS } from '../../../modules/nf-core/bcftools/view/main'
include { CHROMOGRAPH as CHROMOGRAPH_AUTOZYG        } from '../../../modules/nf-core/chromograph/main'
include { RHOCALL_VIZ                               } from '../../../modules/nf-core/rhocall/viz/main'
include { UCSC_WIGTOBIGWIG                          } from '../../../modules/nf-core/ucsc/wigtobigwig/main'

workflow ANNOTATE_RHOCALLVIZ {

    take:
        ch_genome_chrsizes // channel: [mandatory] [ path(sizes) ]
        ch_samples         // channel: [mandatory] [ val(sample_meta) ]
        ch_vcf_tbi         // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ]

    main:
        ch_vcf_tbi
            .combine(ch_samples)
            .map {_meta, vcf, tbi, meta2 -> return [meta2,vcf,tbi]}
            .set { ch_rhocall_viz }

        BCFTOOLS_VIEW_RHOCALL(ch_rhocall_viz, [],[],[])

        BCFTOOLS_VIEW_RHOCALL.out.vcf
            .join(BCFTOOLS_VIEW_RHOCALL.out.tbi)
            .set {ch_roh_in }

        BCFTOOLS_ROH(ch_roh_in, [[],[]], [], [], [], [])

        BCFTOOLS_VIEW_UNCOMPRESS(ch_roh_in,[],[],[])

        BCFTOOLS_VIEW_UNCOMPRESS.out.vcf
                .join(BCFTOOLS_ROH.out.roh)
                .multiMap { meta, vcf, roh ->
                    vcf: [meta, vcf]
                    roh: [meta, roh]
                }
                .set { ch_rhocall_viz_input }

        RHOCALL_VIZ(
            ch_rhocall_viz_input.vcf,
            ch_rhocall_viz_input.roh,
        )

        CHROMOGRAPH_AUTOZYG(RHOCALL_VIZ.out.bed, [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]])

        UCSC_WIGTOBIGWIG(RHOCALL_VIZ.out.wig, ch_genome_chrsizes)

    emit:
        chromograph_autozyg_plots = CHROMOGRAPH_AUTOZYG.out.plots // channel: [ val(meta), path(png) ]
        rhocall_viz_bed           = RHOCALL_VIZ.out.bed           // channel: [ val(meta), path(bed) ]
        rhocall_viz_wig           = RHOCALL_VIZ.out.wig           // channel: [ val(meta), path(wig) ]
        ucsc_wigtobigwig_bw       = UCSC_WIGTOBIGWIG.out.bw       // channel: [ val(meta), path(bw) ]
}
