//
// Annotate MT
//

include { REPLACE_SPACES_IN_VCFINFO                      } from '../../modules/local/replace_spaces_in_vcfinfo'
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
        ch_mt_vcf              // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_mt_tbi              // channel: [mandatory] [ val(meta), path(tbi) ]
        ch_cadd_header         // channel: [mandatory] [ path(txt) ]
        ch_cadd_resources      // channel: [mandatory] [ path(annotation) ]
        ch_genome_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_vcfanno_extra       // channel: [mandatory] [ [path(vcf),path(index)] ]
        ch_vcfanno_lua         // channel: [mandatory] [ path(lua) ]
        ch_vcfanno_resources   // channel: [mandatory] [ [path(vcf1),path(index1),...,path(vcfn),path(indexn)] ]
        ch_vcfanno_toml        // channel: [mandatory] [ path(toml) ]
        val_vep_genome         // string:  [mandatory] GRCh37 or GRCh38
        val_vep_cache_version  // string:  [mandatory] 107
        ch_vep_cache           // channel: [mandatory] [ path(cache) ]
        ch_vep_extra_files     // channel: [mandatory] [ path(files) ]

    main:
        ch_versions     = Channel.empty()
        ch_haplog       = Channel.empty()

        // add prefix to meta
        ch_mt_vcf
            .map { it  ->
                return [it[0]+ [prefix: it[1].simpleName + "_hmtnote"], it[1]]
            }
            .set { ch_hmtnote_in }

        // HMTNOTE ANNOTATE
        HMTNOTE_ANNOTATE(ch_hmtnote_in)
        REPLACE_SPACES_IN_VCFINFO(HMTNOTE_ANNOTATE.out.vcf)
        ZIP_TABIX_HMTNOTE_MT(REPLACE_SPACES_IN_VCFINFO.out.vcf)

        // Vcfanno
        ZIP_TABIX_HMTNOTE_MT.out.gz_tbi
            .combine(ch_vcfanno_extra)
            .map { meta, vcf, tbi, resources -> return [meta + [prefix: meta.prefix + "_vcfanno"], vcf, tbi, resources]}
            .set { ch_in_vcfanno }

        VCFANNO_MT(ch_in_vcfanno, ch_vcfanno_toml, ch_vcfanno_lua, ch_vcfanno_resources)
        ZIP_TABIX_VCFANNO_MT(VCFANNO_MT.out.vcf)

        ch_vcfanno_vcf = ZIP_TABIX_VCFANNO_MT.out.gz_tbi.map{meta, vcf, tbi -> return [meta, vcf]}
        ch_vcfanno_tbi = ZIP_TABIX_VCFANNO_MT.out.gz_tbi.map{meta, vcf, tbi -> return [meta, tbi]}

        // Annotating with CADD
        if (params.cadd_resources != null) {
            ANNOTATE_CADD (
                ZIP_TABIX_VCFANNO_MT.out.gz_tbi,
                ch_cadd_header,
                ch_cadd_resources
            )
            ch_cadd_vcf = ANNOTATE_CADD.out.vcf
            ch_versions = ch_versions.mix(ANNOTATE_CADD.out.versions)
        } else {
            ch_cadd_vcf = Channel.empty()
        }

        // Pick input for vep
        ch_vcfanno_vcf
            .join(ch_cadd_vcf, remainder: true) // If CADD is not run then the third element in this channel will be `null`
            .branch { it  ->                    // If CADD is run, then "it" will be [[meta],selvar.vcf,cadd.vcf], else [[meta],selvar.vcf,null]
                merged: it[2].equals(null)
                    return [it[0]+ [prefix: it[0].prefix + "_vep"], it[1]]
                cadd: !(it[2].equals(null))
                    return [it[0] + [prefix: it[0].prefix + "_cadd_vep"], it[2]]
            }
            .set { ch_for_mix }

        ch_for_mix.merged.mix(ch_for_mix.cadd)
            .tap { ch_haplogrep_in }
            .map { meta, vcf -> return [meta, vcf, []] }
            .set { ch_vep_in }

        // Annotating with ensembl Vep
        ENSEMBLVEP_MT(
            ch_vep_in,
            val_vep_genome,
            "homo_sapiens",
            val_vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta,
            ch_vep_extra_files
        )

        TABIX_TABIX_VEP_MT(ENSEMBLVEP_MT.out.vcf)

        // Running haplogrep3
        if (!(params.skip_tools && params.skip_tools.split(',').contains('haplogrep3'))) {
            HAPLOGREP3_CLASSIFY_MT(ch_haplogrep_in)
            ch_haplog   = HAPLOGREP3_CLASSIFY_MT.out.txt
            ch_versions = ch_versions.mix(HAPLOGREP3_CLASSIFY_MT.out.versions)
        }

        ch_versions = ch_versions.mix(ENSEMBLVEP_MT.out.versions)
        ch_versions = ch_versions.mix(TABIX_TABIX_VEP_MT.out.versions)
        ch_versions = ch_versions.mix(VCFANNO_MT.out.versions)
        ch_versions = ch_versions.mix(HMTNOTE_ANNOTATE.out.versions)
        ch_versions = ch_versions.mix(ZIP_TABIX_VCFANNO_MT.out.versions)
        ch_versions = ch_versions.mix(ZIP_TABIX_HMTNOTE_MT.out.versions)
        ch_versions = ch_versions.mix(REPLACE_SPACES_IN_VCFINFO.out.versions)

    emit:
        haplog    = ch_haplog                   // channel: [ val(meta), path(txt) ]
        vcf_ann   = ENSEMBLVEP_MT.out.vcf       // channel: [ val(meta), path(vcf) ]
        tbi       = TABIX_TABIX_VEP_MT.out.tbi  // channel: [ val(meta), path(tbi) ]
        report    = ENSEMBLVEP_MT.out.report    // channel: [ path(html) ]
        versions  = ch_versions                 // channel: [ path(versions.yml) ]
}
