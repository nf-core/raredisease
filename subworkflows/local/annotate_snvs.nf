//
// A subworkflow to annotate snvs
//

include { VCFANNO                               } from '../../modules/nf-core/vcfanno/main'
include { BCFTOOLS_CONCAT                       } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_ROH                          } from '../../modules/nf-core/bcftools/roh/main'
include { BCFTOOLS_VIEW                         } from '../../modules/nf-core/bcftools/view/main'
include { RHOCALL_ANNOTATE                      } from '../../modules/nf-core/rhocall/annotate/main'
include { ENSEMBLVEP as ENSEMBLVEP_SNV          } from '../../modules/local/ensemblvep/main'
include { TABIX_BGZIPTABIX as ZIP_TABIX_ROHCALL } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as ZIP_TABIX_VCFANNO } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_VEP              } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_BCFTOOLS_CONCAT  } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_BCFTOOLS_VIEW    } from '../../modules/nf-core/tabix/tabix/main'
include { GATK4_SELECTVARIANTS                  } from '../../modules/nf-core/gatk4/selectvariants/main'
include { ANNOTATE_CADD                         } from './annotation/annotate_cadd'

workflow ANNOTATE_SNVS {

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
        ch_fasta              // channel: [mandatory] [ path(fasta) ]
        ch_gnomad_af          // channel: [optional] [ path(tab), path(tbi) ]
        ch_split_intervals    // channel: [mandatory] [ path(intervals) ]
        ch_samples            // channel: [mandatory] [ val(sample_id), val(sex), val(phenotype), val(paternal_id), val(maternal_id), val(case_id) ]

    main:
        ch_versions       = Channel.empty()
        ch_vcf_scatter_in = Channel.empty()
        ch_vep_in         = Channel.empty()

        ch_vcf.map { meta, vcf, idx -> return [vcf, idx] }.set { ch_roh_vcfs }
        ch_samples
            .branch { it ->
                affected: it.phenotype == "2"
                unaffected: it.phenotype == "1"
            }.set { ch_phenotype }
        ch_phenotype.affected.combine(ch_roh_vcfs).set { ch_roh_input }

        BCFTOOLS_ROH (ch_roh_input, ch_gnomad_af, [], [], [], [])

        BCFTOOLS_ROH.out.roh
            .map { meta, roh ->
                new_meta = [:]
                new_meta.id = meta.case_id
                return [new_meta, roh]
            }
            .set { ch_roh_rhocall }

        RHOCALL_ANNOTATE (ch_vcf, ch_roh_rhocall, [])

        ZIP_TABIX_ROHCALL (RHOCALL_ANNOTATE.out.vcf)

        VCFANNO (ZIP_TABIX_ROHCALL.out.gz_tbi, ch_vcfanno_toml, ch_vcfanno_lua, ch_vcfanno_resources)

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
        ANNOTATE_CADD (
            GATK4_SELECTVARIANTS.out.vcf,
            GATK4_SELECTVARIANTS.out.tbi,
            ch_cadd_header,
            ch_cadd_resources
        )

        // If CADD is run, pick CADD output as input for VEP else pass selectvariants output to VEP.
        GATK4_SELECTVARIANTS.out.vcf
            .combine(ANNOTATE_CADD.out.vcf.ifEmpty("null")) // If CADD is not run then this channel will be empty, so assign a default value to allow filtering with branch operator
            .branch { it  ->                                // If CADD is run, then "it" will be [[meta],selvar.vcf,[meta],cadd.vcf], else [[meta],selvar.vcf,null]
                selvar: it[2].equals("null")
                    return [it[0], it[1]]
                cadd: !(it[2].equals("null"))
                    return [it[2], it[3]]
            }
            .set { ch_for_mix }

        ch_vep_in = ch_for_mix.selvar.mix(ch_for_mix.cadd)

        // Annotating with ensembl Vep
        ENSEMBLVEP_SNV(
            ch_vep_in,
            val_vep_genome,
            "homo_sapiens",
            val_vep_cache_version,
            ch_vep_cache,
            ch_fasta,
            []
        )

        ENSEMBLVEP_SNV.out.vcf_gz
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
        ch_versions = ch_versions.mix(ZIP_TABIX_VCFANNO.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)
        ch_versions = ch_versions.mix(TABIX_BCFTOOLS_VIEW.out.versions)
        ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions.first())
        ch_versions = ch_versions.mix(ANNOTATE_CADD.out.versions)
        ch_versions = ch_versions.mix(ENSEMBLVEP_SNV.out.versions.first())
        ch_versions = ch_versions.mix(TABIX_VEP.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
        ch_versions = ch_versions.mix(TABIX_BCFTOOLS_CONCAT.out.versions)

    emit:
        vcf_ann  = ch_vep_ann   // channel: [ val(meta), path(vcf) ]
        tbi      = ch_vep_index // channel: [ val(meta), path(tbi) ]
        versions = ch_versions  // channel: [ path(versions.yml) ]
}
