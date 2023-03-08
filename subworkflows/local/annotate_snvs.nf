//
// A subworkflow to annotate snvs
//

include { VCFANNO                                 } from '../../modules/nf-core/vcfanno/main'
include { BCFTOOLS_CONCAT                         } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_ROH                            } from '../../modules/nf-core/bcftools/roh/main'
include { BCFTOOLS_VIEW                           } from '../../modules/nf-core/bcftools/view/main'
include { RHOCALL_ANNOTATE                        } from '../../modules/nf-core/rhocall/annotate/main'
include { ENSEMBLVEP as ENSEMBLVEP_SNV            } from '../../modules/local/ensemblvep/main'
include { TABIX_BGZIPTABIX as ZIP_TABIX_ROHCALL   } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as ZIP_TABIX_VCFANNO   } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_VEP                } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_BCFTOOLS_CONCAT    } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_BCFTOOLS_VIEW      } from '../../modules/nf-core/tabix/tabix/main'
include { GATK4_SELECTVARIANTS                    } from '../../modules/nf-core/gatk4/selectvariants/main'

workflow ANNOTATE_SNVS {

    take:
        vcf
        vcfanno_resources
        vcfanno_lua
        vcfanno_toml
        vep_genome
        vep_cache_version
        vep_cache
        fasta
        gnomad_af
        split_intervals
        samples

    main:
        ch_versions       = Channel.empty()
        ch_vcf_scatter_in = Channel.empty()
        ch_vep_in         = Channel.empty()

        vcf.map { meta, vcf, idx -> return [vcf, idx] }.set { ch_roh_vcfs }
        samples
            .branch { it ->
                affected: it.phenotype == "2"
                unaffected: it.phenotype == "1"
            }.set { ch_phenotype }
        ch_phenotype.affected.combine(ch_roh_vcfs).set { ch_roh_input }

        BCFTOOLS_ROH (ch_roh_input, gnomad_af, [], [], [], [])

        BCFTOOLS_ROH.out.roh
            .map { meta, roh ->
                new_meta = [:]
                new_meta.id = meta.case_id
                return [new_meta, roh]
            }
            .set { ch_roh_rhocall }

        RHOCALL_ANNOTATE (vcf, ch_roh_rhocall, [])

        ZIP_TABIX_ROHCALL (RHOCALL_ANNOTATE.out.vcf)

        VCFANNO (ZIP_TABIX_ROHCALL.out.gz_tbi, vcfanno_toml, vcfanno_lua, vcfanno_resources)

        ZIP_TABIX_VCFANNO (VCFANNO.out.vcf)

        BCFTOOLS_VIEW(ZIP_TABIX_VCFANNO.out.gz_tbi, [], [], [])  // filter on frequencies

        TABIX_BCFTOOLS_VIEW (BCFTOOLS_VIEW.out.vcf)

        BCFTOOLS_VIEW.out.vcf.join(TABIX_BCFTOOLS_VIEW.out.tbi).collect().set { ch_vcf_scatter_in }

        GATK4_SELECTVARIANTS (ch_vcf_scatter_in.combine(split_intervals)).vcf.set { ch_vep_in }

        ENSEMBLVEP_SNV(
            ch_vep_in,
            vep_genome,
            "homo_sapiens",
            vep_cache_version,
            vep_cache,
            fasta,
            []
        )

        TABIX_VEP (ENSEMBLVEP_SNV.out.vcf_gz)

        ENSEMBLVEP_SNV.out.vcf_gz
            .join(TABIX_VEP.out.tbi)
            .groupTuple()
            .map { meta, vcfs, tbis ->
                def sortedvcfs = vcfs.sort { it.baseName }
                def sortedtbis = tbis.sort { it.baseName }
                return [ meta, sortedvcfs, sortedtbis ]
            }
            .set { ch_vep_ann }

        BCFTOOLS_CONCAT (ch_vep_ann)

        TABIX_BCFTOOLS_CONCAT (BCFTOOLS_CONCAT.out.vcf)

        ch_versions = ch_versions.mix(BCFTOOLS_ROH.out.versions)
        ch_versions = ch_versions.mix(RHOCALL_ANNOTATE.out.versions)
        ch_versions = ch_versions.mix(ZIP_TABIX_ROHCALL.out.versions)
        ch_versions = ch_versions.mix(VCFANNO.out.versions)
        ch_versions = ch_versions.mix(ZIP_TABIX_VCFANNO.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)
        ch_versions = ch_versions.mix(TABIX_BCFTOOLS_VIEW.out.versions)
        ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions.first())
        ch_versions = ch_versions.mix(ENSEMBLVEP_SNV.out.versions.first())
        ch_versions = ch_versions.mix(TABIX_VEP.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)
        ch_versions = ch_versions.mix(TABIX_BCFTOOLS_CONCAT.out.versions)

    emit:
        vcf_ann       = BCFTOOLS_CONCAT.out.vcf
        tbi           = TABIX_BCFTOOLS_CONCAT.out.tbi
        versions      = ch_versions
}
