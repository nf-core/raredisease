//
// Annotate MT
//

include { REPLACE_SPACES_IN_VCFINFO                      } from '../../modules/local/replace_spaces_in_vcfinfo'
include { BCFTOOLS_PLUGINSETGT                           } from '../../modules/nf-core/bcftools/pluginsetgt'
include { TABIX_TABIX as TABIX_TABIX_VEP_MT              } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIPTABIX as ZIP_TABIX_HMTNOTE_MT       } from '../../modules/nf-core/tabix/bgziptabix/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_MT                } from '../../modules/nf-core/ensemblvep/vep/main'
include { HAPLOGREP3_CLASSIFY as HAPLOGREP3_CLASSIFY_MT  } from '../../modules/nf-core/haplogrep3/classify/main'
include { VCFANNO as VCFANNO_MT                          } from '../../modules/nf-core/vcfanno/main'
include { ANNOTATE_CADD                                  } from './annotate_cadd'
include { TABIX_BGZIPTABIX as ZIP_TABIX_VCFANNO_MT       } from '../../modules/nf-core/tabix/bgziptabix/main'
include { HMTNOTE_ANNOTATE                               } from '../../modules/nf-core/hmtnote/annotate/main'

workflow ANNOTATE_MT_SNVS {
    take:
        ch_cadd_header              // channel: [mandatory] [ path(txt) ]
        ch_cadd_resources           // channel: [mandatory] [ path(annotation) ]
        ch_genome_fasta             // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_fai                      // channel: [mandatory] [ path(fai) ]
        ch_mt_vcf                   // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_vcfanno_extra            // channel: [mandatory] [ [path(vcf),path(index)] ]
        ch_vcfanno_lua              // channel: [mandatory] [ path(lua) ]
        ch_vcfanno_resources        // channel: [mandatory] [ [path(vcf1),path(index1),...,path(vcfn),path(indexn)] ]
        ch_vcfanno_toml             // channel: [mandatory] [ path(toml) ]
        ch_vep_cache                // channel: [mandatory] [ path(cache) ]
        ch_vep_extra_files          // channel: [mandatory] [ path(files) ]
        skip_haplogrep3             // boolean
        val_cadd_resources          //  string:  path to cadd resources file
        val_genome                  //  string:  GRCh37 or GRCh38
        val_homoplasmy_af_threshold //   float: 0-1
        val_vep_cache_version       //  string:  vep version ex: 107

    main:
        ch_versions     = channel.empty()
        ch_haplog       = channel.empty()

        // add prefix to meta
        ch_mt_vcf
            .map { meta, vcf  ->
                return [meta+ [prefix: vcf.simpleName + "_hmtnote"], vcf]
            }
            .set { ch_hmtnote_in }

        // HMTNOTE ANNOTATE
        HMTNOTE_ANNOTATE(ch_hmtnote_in)
        REPLACE_SPACES_IN_VCFINFO(HMTNOTE_ANNOTATE.out.vcf)
        ZIP_TABIX_HMTNOTE_MT(REPLACE_SPACES_IN_VCFINFO.out.vcf)

        // Vcfanno
        ZIP_TABIX_HMTNOTE_MT.out.gz_index
            .combine(ch_vcfanno_extra)
            .map { meta, vcf, tbi, resources -> return [meta + [prefix: meta.prefix + "_vcfanno"], vcf, tbi, resources]}
            .set { ch_in_vcfanno }

        VCFANNO_MT(ch_in_vcfanno, ch_vcfanno_toml, ch_vcfanno_lua, ch_vcfanno_resources)
        ZIP_TABIX_VCFANNO_MT(VCFANNO_MT.out.vcf)

        ch_vcfanno_vcf = ZIP_TABIX_VCFANNO_MT.out.gz_index.map{meta, vcf, _tbi -> return [meta, vcf]}

        // Annotating with CADD
        if (!val_cadd_resources.equals(null)) {
            ANNOTATE_CADD (
                ch_cadd_resources,
                ch_fai,
                ch_cadd_header,
                ZIP_TABIX_VCFANNO_MT.out.gz_index,
                val_genome
            )
            ch_cadd_vcf = ANNOTATE_CADD.out.vcf
            ch_versions = ch_versions.mix(ANNOTATE_CADD.out.versions)
        } else {
            ch_cadd_vcf = channel.empty()
        }

        ch_vcfanno_vcf
            .join(ch_cadd_vcf, remainder: true)
            .branch { meta, vcfanno, cadd  ->
                vcfanno: cadd.equals(null)
                    return [meta+ [prefix: meta.prefix + "_vep"], vcfanno]
                cadd: !(cadd.equals(null))
                    return [meta + [prefix: meta.prefix + "_cadd_vep"], cadd]
            }
            .set { ch_annotated_vcfs }

        ch_annotated_vcfs.vcfanno.mix(ch_annotated_vcfs.cadd)
            .tap { ch_haplogrep_in }
            .map { meta, vcf -> return [meta, vcf, []] }
            .set { ch_vep_in }

        // Annotating with ensembl Vep
        ENSEMBLVEP_MT(
            ch_vep_in,
            val_genome,
            "homo_sapiens",
            val_vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta,
            ch_vep_extra_files
        )

        TABIX_TABIX_VEP_MT(ENSEMBLVEP_MT.out.vcf)

        ch_vcf = ENSEMBLVEP_MT.out.vcf
        ch_tbi = TABIX_TABIX_VEP_MT.out.index

        // Running haplogrep3
        if (!skip_haplogrep3) {
            HAPLOGREP3_CLASSIFY_MT(ch_haplogrep_in)
            ch_haplog   = HAPLOGREP3_CLASSIFY_MT.out.txt
        }

        if (val_homoplasmy_af_threshold<1) {
            BCFTOOLS_PLUGINSETGT (
                ENSEMBLVEP_MT.out.vcf.map { meta, vcf -> return [meta, vcf, []] },
                channel.value('q'),
                channel.value("c:'1/1'"),
                [],
                []
            )
            ch_vcf = BCFTOOLS_PLUGINSETGT.out.vcf
            ch_tbi = BCFTOOLS_PLUGINSETGT.out.tbi
        }

        ch_versions = ch_versions.mix(HMTNOTE_ANNOTATE.out.versions)
        ch_versions = ch_versions.mix(REPLACE_SPACES_IN_VCFINFO.out.versions)

    emit:
        haplog    = ch_haplog                // channel: [ val(meta), path(txt) ]
        report    = ENSEMBLVEP_MT.out.report // channel: [ path(html) ]
        tbi       = ch_tbi                   // channel: [ val(meta), path(tbi) ]
        vcf_ann   = ch_vcf                   // channel: [ val(meta), path(vcf) ]
        versions  = ch_versions              // channel: [ path(versions.yml) ]
}
