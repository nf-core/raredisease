//
// A subworkflow to annotate snvs
//

include { VCFANNO                        } from '../../modules/nf-core/modules/vcfanno/main'
include { BCFTOOLS_ROH                   } from '../../modules/nf-core/modules/bcftools/roh/main'
include { RHOCALL_ANNOTATE               } from '../../modules/nf-core/modules/rhocall/annotate/main'
include { ENSEMBLVEP as ENSEMBLVEP_SNV   } from '../../modules/local/ensemblvep/main'
include { TABIX_TABIX as TABIX_SNV_ANNO  } from '../../modules/nf-core/modules/tabix/tabix/main'


workflow ANNOTATE_SNVS {

    take:
        vcf
        vcfanno_resource_dir
        vcfanno_toml
        vep_genome
        vep_cache_version
        vep_cache
        fasta
        gnomad_af
        samples

    main:
        ch_versions = Channel.empty()
        ch_toml     = file(vcfanno_toml)

        //
        // annotate vcfanno
        //

        VCFANNO (vcf, ch_toml, vcfanno_resource_dir)
        ch_versions = ch_versions.mix(VCFANNO.out.versions)

        //
        // annotate rhocall
        //
        vcf.map { meta, vcf, idx ->
                return [ vcf, idx ]
            }
            .set { ch_roh_vcfs}

        samples
            .branch { it ->
                affected: it.phenotype == "2"
                unaffected: it.phenotype == "1"
            }
            .set {ch_phenotype}

        ch_phenotype.affected
            .combine(ch_roh_vcfs)
            .set { ch_roh_input }

        BCFTOOLS_ROH (ch_roh_input, gnomad_af, [], [], [], [])

        BCFTOOLS_ROH.out.roh
            .map { meta, roh ->
                new_meta = [:]
                new_meta.id = meta.case_id
                return [new_meta, roh]
            }
            .set { ch_roh_rhocall}

        RHOCALL_ANNOTATE (vcf, ch_roh_rhocall, [])

        TABIX_SNV_ANNO (RHOCALL_ANNOTATE.out.vcf)
        ch_versions = ch_versions.mix(TABIX_SNV_ANNO.out.versions)

        RHOCALL_ANNOTATE.out
            .vcf
            .join(TABIX_SNV_ANNO.out.tbi)
            .set { ch_vep_in }

        ENSEMBLVEP_SNV(ch_vep_in,
            vep_genome,
            "homo_sapiens",
            vep_cache_version,
            vep_cache
            )
        ch_versions = ch_versions.mix(ENSEMBLVEP_SNV.out.versions)

    emit:
        vcf_ann                = ENSEMBLVEP_SNV.out.vcf
        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
