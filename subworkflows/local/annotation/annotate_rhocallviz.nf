//
// A subworkflow to plot binned zygosity and RHO-regions.
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_RHOCALL    } from '../../../modules/nf-core/bcftools/view/main'
include { TABIX_TABIX                               } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_ROH                              } from '../../../modules/nf-core/bcftools/roh/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_UNCOMPRESS } from '../../../modules/nf-core/bcftools/view/main'
include { RHOCALL_VIZ                               } from '../../../modules/nf-core/rhocall/viz/main'
include { UCSC_WIGTOBIGWIG                          } from '../../../modules/nf-core/ucsc/wigtobigwig/main'
include { CHROMOGRAPH as CHROMOGRAPH_AUTOZYG        } from '../../../modules/nf-core/chromograph/main'

workflow ANNOTATE_RHOCALLVIZ {

    take:
        ch_vcf_tbi         // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ]
        ch_samples         // channel: [mandatory] [ val(sample_meta) ]
        ch_genome_chrsizes // channel: [mandatory] [ path(sizes) ]

    main:
        ch_versions       = Channel.empty()

        ch_vcf_tbi
            .combine(ch_samples)
            .map {meta, vcf, tbi, meta2 -> return [meta2,vcf,tbi]}
            .set { ch_rhocall_viz }

        BCFTOOLS_VIEW_RHOCALL(ch_rhocall_viz, [],[],[])

        TABIX_TABIX(BCFTOOLS_VIEW_RHOCALL.out.vcf)

        BCFTOOLS_VIEW_RHOCALL.out.vcf
            .join(TABIX_TABIX.out.tbi)
            .set {ch_roh_in }

        BCFTOOLS_ROH(ch_roh_in, [[],[]], [], [], [], [])

        BCFTOOLS_VIEW_UNCOMPRESS(ch_roh_in,[],[],[])

        RHOCALL_VIZ(BCFTOOLS_VIEW_UNCOMPRESS.out.vcf, BCFTOOLS_ROH.out.roh)

        CHROMOGRAPH_AUTOZYG(RHOCALL_VIZ.out.bed, [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]])

        UCSC_WIGTOBIGWIG(RHOCALL_VIZ.out.wig, ch_genome_chrsizes)

        ch_versions = ch_versions.mix(BCFTOOLS_VIEW_RHOCALL.out.versions.first())
        ch_versions = ch_versions.mix(CHROMOGRAPH_AUTOZYG.out.versions.first())
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_ROH.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW_UNCOMPRESS.out.versions.first())
        ch_versions = ch_versions.mix(RHOCALL_VIZ.out.versions.first())
        ch_versions = ch_versions.mix(UCSC_WIGTOBIGWIG.out.versions.first())

    emit:
        versions = ch_versions  // channel: [ path(versions.yml) ]
}
