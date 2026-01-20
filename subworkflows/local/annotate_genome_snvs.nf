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
include { ANNOTATE_CADD                         } from './annotate_cadd'
include { ANNOTATE_RHOCALLVIZ                   } from './annotate_rhocallviz'

workflow ANNOTATE_GENOME_SNVS {

    take:
        ch_cadd_header        // channel: [mandatory] [ path(txt) ]
        ch_cadd_resources     // channel: [mandatory] [ path(annotation) ]
        ch_genome_chrsizes    // channel: [mandatory] [ path(sizes) ]
        ch_genome_fai         // channel: [mandatory] [ path(fai) ]
        ch_genome_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_gnomad_af          // channel: [optional] [ path(tab), path(tbi) ]
        ch_samples            // channel: [mandatory] [ val(sample_meta) ]
        ch_split_intervals    // channel: [mandatory] [ path(intervals) ]
        ch_vcf                // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ]
        ch_vcfanno_extra      // channel: [mandatory] [ [path(vcf),path(index)] ]
        ch_vcfanno_lua        // channel: [mandatory] [ path(lua) ]
        ch_vcfanno_resources  // channel: [mandatory] [ [path(vcf1),path(index1),...,path(vcfn),path(indexn)] ]
        ch_vcfanno_toml       // channel: [mandatory] [ path(toml) ]
        ch_vep_cache          // channel: [mandatory] [ path(cache) ]
        ch_vep_extra_files    // channel: [mandatory] [ path(files) ]
        val_cadd_resources    // string: path to cadd resources file
        val_genome            // string: GRCh37 or GRCh38
        val_vep_cache_version // string:  vep version ex: 107

    main:
        ch_cadd_vcf       = channel.empty()
        ch_versions       = channel.empty()
        ch_vcf_scatter_in = channel.empty()
        ch_vep_in         = channel.empty()

        ch_vcf
            .filter { meta, _vcf, _tbi ->
                meta.probands.size() > 0
            }
            .set { ch_roh_in }

        BCFTOOLS_ROH (ch_roh_in, ch_gnomad_af, [], [], [], [])

        RHOCALL_ANNOTATE (ch_vcf, BCFTOOLS_ROH.out.roh, [])

        ZIP_TABIX_ROHCALL (RHOCALL_ANNOTATE.out.vcf)

        ch_vcf
            .join(ZIP_TABIX_ROHCALL.out.gz_tbi, remainder: true)
            .combine(ch_split_intervals)
            .map { it ->
                    def meta = it[0]
                    def vcf  = it[1]
                    def tbi  = it[2]

                def hasRohCall = (it.size() == 6)

                if (hasRohCall) {
                    def rohcall      = it[3]
                    def rohcallindex = it[4]
                    def interval     = it[5]
                    return [
                        meta + [prefix: meta.id + "_rhocall", scatterid: interval.baseName],
                        rohcall,
                        rohcallindex,
                        interval
                    ]
                } else {
                    def interval = it[4]
                    return [
                        meta + [prefix: meta.id, scatterid: interval.baseName],
                        vcf,
                        tbi,
                        interval
                    ]
                }
            }
            .set { ch_vcf_scatter_in }

        GATK4_SELECTVARIANTS (ch_vcf_scatter_in)

        GATK4_SELECTVARIANTS.out.vcf
            .join(GATK4_SELECTVARIANTS.out.tbi)
            .combine(ch_vcfanno_extra)
            .set { ch_vcfanno_in }

        VCFANNO (ch_vcfanno_in, ch_vcfanno_toml, ch_vcfanno_lua, ch_vcfanno_resources)

        ZIP_TABIX_VCFANNO (VCFANNO.out.vcf)

        BCFTOOLS_VIEW(ZIP_TABIX_VCFANNO.out.gz_tbi, [], [], [])  // filter on frequencies

        // Annotating with CADD
        if (!val_cadd_resources.equals(null)) {
            TABIX_BCFTOOLS_VIEW (BCFTOOLS_VIEW.out.vcf)

            BCFTOOLS_VIEW.out.vcf
                .join(TABIX_BCFTOOLS_VIEW.out.tbi, failOnMismatch:true, failOnDuplicate:true)
                .set { ch_cadd_in }

            ANNOTATE_CADD (
                ch_cadd_resources,
                ch_genome_fai,
                ch_cadd_header,
                ch_cadd_in,
                val_genome
            )
            ch_cadd_vcf = ANNOTATE_CADD.out.vcf
            ch_versions = ch_versions.mix(ANNOTATE_CADD.out.versions)
            ch_versions = ch_versions.mix(TABIX_BCFTOOLS_VIEW.out.versions)
        }

        BCFTOOLS_VIEW.out.vcf
            .join(ch_cadd_vcf, remainder: true)
            .branch { meta, selectvariants, cadd  ->
                selvar: cadd.equals(null)
                    return [meta + [prefix: meta.prefix + "_filter"], selectvariants]
                cadd: !(cadd.equals(null))
                    return [meta + [prefix: meta.prefix + "_filter_cadd"], cadd]
            }
            .set { ch_annotated_vcfs }

        ch_annotated_vcfs.selvar.mix(ch_annotated_vcfs.cadd)
            .map { meta, vcf -> return [meta, vcf, []] }
            .set { ch_vep_in }

        // Annotating with ensembl Vep
        ENSEMBLVEP_SNV(
            ch_vep_in,
            val_genome,
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
                def sortedvcfs = vcfs.sort { vcf -> vcf.baseName }
                def sortedtbis = tbis.sort { tbi -> tbi.baseName }
                return [ meta, sortedvcfs, sortedtbis ]
            }
            .set { ch_concat_in }

        BCFTOOLS_CONCAT (ch_concat_in)

        BCFTOOLS_CONCAT.out.vcf
            .flatMap { meta, vcf ->
                meta.upd_children.collect { upd_sample ->
                    def new_meta = meta + [upd_child: upd_sample, prefix: meta.prefix + "_vcfanno"]
                    [new_meta, vcf]
                }
            }
            .set { ch_upd_in }

        UPD_SITES(ch_upd_in)
        UPD_REGIONS(ch_upd_in)
        CHROMOGRAPH_SITES([[],[]], [[],[]], [[],[]], [[],[]], [[],[]], [[],[]], UPD_SITES.out.bed)
        CHROMOGRAPH_REGIONS([[],[]], [[],[]], [[],[]], [[],[]], [[],[]], UPD_REGIONS.out.bed, [[],[]])


        BCFTOOLS_CONCAT.out.vcf
            .map { meta, vcf -> [meta - meta.subMap('prefix'), vcf] }
            .set { ch_concat_out }

        TABIX_BCFTOOLS_CONCAT (ch_concat_out)

        ch_vep_ann       = ch_concat_out
        ch_vep_index     = TABIX_BCFTOOLS_CONCAT.out.tbi
        ch_vep_ann_index = ch_concat_out.join(TABIX_BCFTOOLS_CONCAT.out.tbi)
        //rhocall_viz
        ANNOTATE_RHOCALLVIZ(ch_genome_chrsizes, ch_samples, ch_vep_ann_index )

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
        ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions.first())
        ch_versions = ch_versions.mix(ENSEMBLVEP_SNV.out.versions.first())
        ch_versions = ch_versions.mix(TABIX_VEP.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
        ch_versions = ch_versions.mix(TABIX_BCFTOOLS_CONCAT.out.versions)
        ch_versions = ch_versions.mix(ANNOTATE_RHOCALLVIZ.out.versions)

    emit:
        tbi      = ch_vep_index // channel: [ val(meta), path(tbi) ]
        vcf_ann  = ch_vep_ann   // channel: [ val(meta), path(vcf) ]
        versions = ch_versions  // channel: [ path(versions.yml) ]
}
