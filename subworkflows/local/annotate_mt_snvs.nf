//
// Annotate MT
//

include { TABIX_TABIX as TABIX_TABIX_MT                  } from '../../modules/nf-core/tabix/tabix/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_MT                } from '../../modules/nf-core/ensemblvep/vep/main'
include { HAPLOGREP2_CLASSIFY as HAPLOGREP2_CLASSIFY_MT  } from '../../modules/nf-core/haplogrep2/classify/main'
include { VCFANNO as VCFANNO_MT                          } from '../../modules/nf-core/vcfanno/main'
include { ANNOTATE_CADD                                  } from './annotation/annotate_cadd'
include { TABIX_BGZIPTABIX as ZIP_TABIX_HMTNOTE          } from '../../modules/nf-core/tabix/bgziptabix/main'
include { HMTNOTE_ANNOTATE                               } from '../../modules/nf-core/hmtnote/annotate/main'

workflow ANNOTATE_MT_SNVS {
    take:
        ch_mt_vcf              // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_mt_tbi              // channel: [mandatory] [ val(meta), path(tbi) ]
        ch_cadd_header         // channel: [mandatory] [ path(txt) ]
        ch_cadd_resources      // channel: [mandatory] [ path(annotation) ]
        ch_genome_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_vcfanno_resources   // channel: [mandatory] [ path(resources) ]
        ch_vcfanno_toml        // channel: [mandatory] [ path(toml) ]
        val_vep_genome         // string:  [mandatory] GRCh37 or GRCh38
        val_vep_cache_version  // string:  [mandatory] 107
        ch_vep_cache           // channel: [mandatory] [ path(cache) ]
        ch_vep_cache           // channel: [mandatory] [ path(cache) ]
        ch_vep_extra_files     // channel: [mandatory] [ path(files) ]

    main:
        ch_cadd_vcf = Channel.empty()
        ch_versions = Channel.empty()

        // Annotating with CADD
        if (params.cadd_resources != null) {
            ANNOTATE_CADD (
                ch_mt_vcf,
                ch_mt_tbi,
                ch_cadd_header,
                ch_cadd_resources
            )
            ch_cadd_vcf = ANNOTATE_CADD.out.vcf
            ch_versions = ch_versions.mix(ANNOTATE_CADD.out.versions)
        }

        // Pick input for vep
        ch_mt_vcf
            .join(ch_cadd_vcf, remainder: true) // If CADD is not run then the third element in this channel will be `null`
            .branch { it  ->                    // If CADD is run, then "it" will be [[meta],selvar.vcf,cadd.vcf], else [[meta],selvar.vcf,null]
                merged: it[2].equals(null)
                    return [it[0], it[1]]
                cadd: !(it[2].equals(null))
                    return [it[0], it[2]]
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

        // Running vcfanno
        TABIX_TABIX_MT(ENSEMBLVEP_MT.out.vcf)
        ENSEMBLVEP_MT.out.vcf
            .join(TABIX_TABIX_MT.out.tbi, failOnMismatch:true, failOnDuplicate:true)
            .map { meta, vcf, tbi -> return [meta, vcf, tbi, []]}
            .set { ch_in_vcfanno }

        VCFANNO_MT(ch_in_vcfanno, ch_vcfanno_toml, [], ch_vcfanno_resources)

        // HMTNOTE ANNOTATE
        HMTNOTE_ANNOTATE(VCFANNO_MT.out.vcf)
        HMTNOTE_ANNOTATE.out.vcf.map{meta, vcf ->
            return [meta, WorkflowRaredisease.replaceSpacesInInfoColumn(vcf, vcf.parent.toString(), vcf.baseName)]
            }
            .set { ch_hmtnote_reformatted }
        ZIP_TABIX_HMTNOTE(ch_hmtnote_reformatted)

        // Prepare output
        ch_vcf_out = ZIP_TABIX_HMTNOTE.out.gz_tbi.map{meta, vcf, tbi -> return [meta, vcf] }
        ch_tbi_out = ZIP_TABIX_HMTNOTE.out.gz_tbi.map{meta, vcf, tbi -> return [meta, tbi] }

        // Running haplogrep2
        HAPLOGREP2_CLASSIFY_MT(ch_haplogrep_in, "vcf.gz")

        ch_versions = ch_versions.mix(ENSEMBLVEP_MT.out.versions)
        ch_versions = ch_versions.mix(TABIX_TABIX_MT.out.versions)
        ch_versions = ch_versions.mix(VCFANNO_MT.out.versions)
        ch_versions = ch_versions.mix(HMTNOTE_ANNOTATE.out.versions)
        ch_versions = ch_versions.mix(HAPLOGREP2_CLASSIFY_MT.out.versions)

    emit:
        haplog    = HAPLOGREP2_CLASSIFY_MT.out.txt // channel: [ val(meta), path(txt) ]
        vcf_ann   = ch_vcf_out                     // channel: [ val(meta), path(vcf) ]
        tbi       = ch_tbi_out                     // channel: [ val(meta), path(tbi) ]
        report    = ENSEMBLVEP_MT.out.report       // channel: [ path(html) ]
        versions  = ch_versions                    // channel: [ path(versions.yml) ]
}
