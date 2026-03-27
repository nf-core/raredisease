//
// Annotate MT
//


include { BCFTOOLS_PLUGINSETGT                           } from '../../../modules/nf-core/bcftools/pluginsetgt'
include { TABIX_TABIX as TABIX_TABIX_VEP_MT              } from '../../../modules/nf-core/tabix/tabix/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_MT                } from '../../../modules/nf-core/ensemblvep/vep/main'
include { HAPLOGREP3_CLASSIFY as HAPLOGREP3_CLASSIFY_MT  } from '../../../modules/nf-core/haplogrep3/classify/main'
include { VCFANNO as VCFANNO_MT                          } from '../../../modules/nf-core/vcfanno/main'
include { ANNOTATE_CADD                                  } from '../annotate_cadd'

workflow ANNOTATE_MT_SNVS {
    take:
        ch_cadd_header              // channel: [mandatory] [ path(txt) ]
        ch_cadd_resources           // channel: [mandatory] [ path(annotation) ]
        ch_genome_fasta             // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_fai                      // channel: [mandatory] [ path(fai) ]
        ch_mt_vcf_tbi               // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ]
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
        ch_haplog       = channel.empty()

        // Vcfanno
        ch_mt_vcf_tbi
            .combine(ch_vcfanno_extra)
            .map { meta, vcf, tbi, resources -> return [meta + [prefix: vcf.simpleName + "_vcfanno"], vcf, tbi, resources]}
            .set { ch_in_vcfanno }

        VCFANNO_MT(ch_in_vcfanno, ch_vcfanno_toml, ch_vcfanno_lua, ch_vcfanno_resources)

        // Annotating with CADD
        if (!val_cadd_resources.equals(null)) {
            ANNOTATE_CADD (
                ch_cadd_resources,
                ch_fai,
                ch_cadd_header,
                VCFANNO_MT.out.tbi,
                val_genome
            )
            ch_cadd_vcf = ANNOTATE_CADD.out.vcf
        } else {
            ch_cadd_vcf = channel.empty()
        }

        VCFANNO_MT.out.vcf
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

        ch_vcf = ENSEMBLVEP_MT.out.vcf
        ch_tbi = ENSEMBLVEP_MT.out.tbi

        // Running haplogrep3
        ch_haplog_publish = channel.empty()
        if (!skip_haplogrep3) {
            HAPLOGREP3_CLASSIFY_MT(ch_haplogrep_in)
            ch_haplog         = HAPLOGREP3_CLASSIFY_MT.out.txt
            ch_haplog_publish = HAPLOGREP3_CLASSIFY_MT.out.txt
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

        ch_publish = ENSEMBLVEP_MT.out.vcf
            .mix(ENSEMBLVEP_MT.out.tbi)
            .mix(ch_haplog_publish)
            .map { meta, value -> ['annotate_snv/mitochondria/', [meta, value]] }

    emit:
        haplog    = ch_haplog                // channel: [ val(meta), path(txt) ]
        publish   = ch_publish               // channel: [ val(destination), val(value) ]
        report    = ENSEMBLVEP_MT.out.report // channel: [ path(html) ]
        tbi       = ch_tbi                   // channel: [ val(meta), path(tbi) ]
        vcf_ann   = ch_vcf                   // channel: [ val(meta), path(vcf) ]
}
