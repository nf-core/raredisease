//
// A subworkflow to annotate snvs in the genome
//

include { VCFANNO                               } from '../../modules/nf-core/vcfanno/main'
include { BCFTOOLS_CONCAT                       } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_ROH                          } from '../../modules/nf-core/bcftools/roh/main'
include { BCFTOOLS_VIEW                         } from '../../modules/nf-core/bcftools/view/main'
include { RHOCALL_ANNOTATE                      } from '../../modules/nf-core/rhocall/annotate/main'
include { UPD as UPD_SITES                      } from '../../modules/nf-core/upd/main'
include { UPD as UPD_REGIONS                    } from '../../modules/nf-core/upd/main'
include { CHROMOGRAPH as CHROMOGRAPH_SITES      } from '../../modules/nf-core/chromograph/main'
include { CHROMOGRAPH as CHROMOGRAPH_REGIONS    } from '../../modules/nf-core/chromograph/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_SNV      } from '../../modules/nf-core/ensemblvep/vep/main'
include { TABIX_BGZIPTABIX as ZIP_TABIX_ROHCALL } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as ZIP_TABIX_VCFANNO } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_VEP              } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_BCFTOOLS_CONCAT  } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_BCFTOOLS_VIEW    } from '../../modules/nf-core/tabix/tabix/main'
include { GATK4_SELECTVARIANTS                  } from '../../modules/nf-core/gatk4/selectvariants/main'
include { ANNOTATE_CADD                         } from './annotation/annotate_cadd'

workflow ANNOTATE_GENOME_SNVS {

    take:
        ch_vcf                // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ]
        analysis_type         // string: [mandatory] 'wgs' or 'wes'
        ch_cadd_header        // channel: [mandatory] [ path(txt) ]
        ch_cadd_resources     // channel: [mandatory] [ path(annotation) ]
        ch_vcfanno_resources  // channel: [mandatory] [ path(resources) ]
        ch_vcfanno_lua        // channel: [mandatory] [ path(lua) ]
        ch_vcfanno_toml       // channel: [mandatory] [ path(toml) ]
        val_vep_genome        // string: [mandatory] GRCh37 or GRCh38
        val_vep_cache_version // string: [mandatory] default: 107
        ch_vep_cache          // channel: [mandatory] [ path(cache) ]
        ch_genome_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_gnomad_af          // channel: [optional] [ path(tab), path(tbi) ]
        ch_split_intervals    // channel: [mandatory] [ path(intervals) ]
        ch_vep_extra_files    // channel: [mandatory] [ path(files) ]

    main:
        ch_cadd_vcf       = Channel.empty()
        ch_versions       = Channel.empty()
        ch_vcf_scatter_in = Channel.empty()
        ch_vep_in         = Channel.empty()

        BCFTOOLS_ROH (ch_vcf, ch_gnomad_af, [], [], [], [])

        RHOCALL_ANNOTATE (ch_vcf, BCFTOOLS_ROH.out.roh, [])

        ZIP_TABIX_ROHCALL (RHOCALL_ANNOTATE.out.vcf)

        ZIP_TABIX_ROHCALL.out.gz_tbi
            .map { meta, vcf, tbi -> return [meta, vcf, tbi, []]}
            .set { ch_vcf_in }

        VCFANNO (ch_vcf_in, ch_vcfanno_toml, ch_vcfanno_lua, ch_vcfanno_resources)

        VCFANNO.out.vcf
            .map {meta, vcf ->
                def splitchannels = []
                for (int i=0; i< meta.upd_children.size(); i++) {
                    upd_sample = meta.upd_children[i]
                    new_meta = meta + [upd_child:upd_sample]
                    splitchannels.add([new_meta,vcf])
                    }
                return splitchannels
            }
            .flatten()
            .buffer (size: 2)
            .set { ch_upd_in }

        UPD_SITES(ch_upd_in)
        UPD_REGIONS(ch_upd_in)
        CHROMOGRAPH_SITES([[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], UPD_SITES.out.bed)
        CHROMOGRAPH_REGIONS([[],[]], [[],[]], [[],[]], [[],[]], [[],[]], UPD_REGIONS.out.bed, [[],[]])

        ZIP_TABIX_VCFANNO (VCFANNO.out.vcf)

        BCFTOOLS_VIEW(ZIP_TABIX_VCFANNO.out.gz_tbi, [], [], [])  // filter on frequencies

        TABIX_BCFTOOLS_VIEW (BCFTOOLS_VIEW.out.vcf)

        BCFTOOLS_VIEW.out.vcf
            .join(TABIX_BCFTOOLS_VIEW.out.tbi, failOnMismatch:true, failOnDuplicate:true)
            .collect()
            .combine(ch_split_intervals)
            .map {
                meta, vcf, tbi, interval ->
                return [meta + [scatterid:interval.baseName], vcf, tbi, interval]
            }
            .set { ch_vcf_scatter_in }

        GATK4_SELECTVARIANTS (ch_vcf_scatter_in)

        // Annotating with CADD
        if (params.cadd_resources != null) {
            ANNOTATE_CADD (
                GATK4_SELECTVARIANTS.out.vcf,
                GATK4_SELECTVARIANTS.out.tbi,
                ch_cadd_header,
                ch_cadd_resources
            )
            ch_cadd_vcf = ANNOTATE_CADD.out.vcf
            ch_versions = ch_versions.mix(ANNOTATE_CADD.out.versions)
        }

        // If CADD is run, pick CADD output as input for VEP else pass selectvariants output to VEP.
        GATK4_SELECTVARIANTS.out.vcf
            .join(ch_cadd_vcf, remainder: true) // If CADD is not run then the third element in this channel will be `null`
            .branch { it  ->                              // If CADD is run, then "it" will be [[meta],selvar.vcf,cadd.vcf], else [[meta],selvar.vcf,null]
                selvar: it[2].equals(null)
                    return [it[0], it[1]]
                cadd: !(it[2].equals(null))
                    return [it[0], it[2]]
            }
            .set { ch_for_mix }

        ch_for_mix.selvar.mix(ch_for_mix.cadd)
            .map { meta, vcf -> return [meta, vcf, []] }
            .set { ch_vep_in }


        // Annotating with ensembl Vep
        ENSEMBLVEP_SNV(
            ch_vep_in,
            val_vep_genome,
            "homo_sapiens",
            val_vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta,
            ch_vep_extra_files
        )

        ENSEMBLVEP_SNV.out.vcf
            .map { meta, vcf -> [meta - meta.subMap('scatterid'), vcf] }
            .set { ch_vep_out }

        TABIX_VEP (ch_vep_out)

        ch_vep_out
            .join(TABIX_VEP.out.tbi, failOnMismatch:true)
            .groupTuple()
            .map { meta, vcfs, tbis ->
                def sortedvcfs = vcfs.sort { it.baseName }
                def sortedtbis = tbis.sort { it.baseName }
                return [ meta, sortedvcfs, sortedtbis ]
            }
            .set { ch_concat_in }

        BCFTOOLS_CONCAT (ch_concat_in)

        TABIX_BCFTOOLS_CONCAT (BCFTOOLS_CONCAT.out.vcf)

        ch_vep_ann   = BCFTOOLS_CONCAT.out.vcf
        ch_vep_index = TABIX_BCFTOOLS_CONCAT.out.tbi

        ch_versions = ch_versions.mix(BCFTOOLS_ROH.out.versions)
        ch_versions = ch_versions.mix(RHOCALL_ANNOTATE.out.versions)
        ch_versions = ch_versions.mix(ZIP_TABIX_ROHCALL.out.versions)
        ch_versions = ch_versions.mix(VCFANNO.out.versions)
        ch_versions = ch_versions.mix(UPD_SITES.out.versions)
        ch_versions = ch_versions.mix(UPD_REGIONS.out.versions)
        ch_versions = ch_versions.mix(CHROMOGRAPH_SITES.out.versions)
        ch_versions = ch_versions.mix(CHROMOGRAPH_REGIONS.out.versions)
        ch_versions = ch_versions.mix(ZIP_TABIX_VCFANNO.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)
        ch_versions = ch_versions.mix(TABIX_BCFTOOLS_VIEW.out.versions)
        ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions.first())
        ch_versions = ch_versions.mix(ENSEMBLVEP_SNV.out.versions.first())
        ch_versions = ch_versions.mix(TABIX_VEP.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
        ch_versions = ch_versions.mix(TABIX_BCFTOOLS_CONCAT.out.versions)

    emit:
        vcf_ann  = ch_vep_ann   // channel: [ val(meta), path(vcf) ]
        tbi      = ch_vep_index // channel: [ val(meta), path(tbi) ]
        versions = ch_versions  // channel: [ path(versions.yml) ]
}
