//
// A subworkflow to annotate snvs
//

include { VCFANNO                             } from '../../modules/nf-core/vcfanno/main'
include { BCFTOOLS_ROH                        } from '../../modules/nf-core/bcftools/roh/main'
include { BCFTOOLS_VIEW                       } from '../../modules/nf-core/bcftools/view/main'
include { RHOCALL_ANNOTATE                    } from '../../modules/nf-core/rhocall/annotate/main'
include { ENSEMBLVEP as ENSEMBLVEP_SNV        } from '../../modules/local/ensemblvep/main'
include { TABIX_BGZIPTABIX as TABIX_ROHCALL   } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_VCFANNO   } from '../../modules/nf-core/tabix/bgziptabix/main'


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
        ch_toml     = Channel.fromPath(vcfanno_toml)

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

        TABIX_ROHCALL (RHOCALL_ANNOTATE.out.vcf)
        ch_versions = ch_versions.mix(TABIX_ROHCALL.out.versions)

        //
        // annotate vcfanno
        //
        VCFANNO (TABIX_ROHCALL.out.gz_tbi, ch_toml, vcfanno_resource_dir)
        ch_versions = ch_versions.mix(VCFANNO.out.versions)

        TABIX_VCFANNO (VCFANNO.out.vcf)
        ch_versions = ch_versions.mix(TABIX_VCFANNO.out.versions)

        BCFTOOLS_VIEW(TABIX_VCFANNO.out.gz_tbi,[],[],[])
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

        //
        // annotate vep
        //
        ENSEMBLVEP_SNV(BCFTOOLS_VIEW.out.vcf,
            vep_genome,
            "homo_sapiens",
            vep_cache_version,
            vep_cache,
            fasta,
            []
            )
        ch_versions = ch_versions.mix(ENSEMBLVEP_SNV.out.versions)

    emit:
        vcf_ann       = ENSEMBLVEP_SNV.out.vcf
        versions      = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
