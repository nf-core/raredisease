//
// A subworkflow to annotate snvs
//

include { VCFANNO      } from '../../modules/nf-core/modules/vcfanno/main'
include { BCFTOOLS_ROH } from '../../modules/nf-core/modules/bcftools/roh/main'


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
        vcf.map { meta, vcf, idx ->
                    return [meta, []]
            }
            .set { ch_placeholder }

        VCFANNO (vcf, ch_placeholder, ch_toml, vcfanno_resource_dir)
        ch_versions = ch_versions.mix(VCFANNO.out.versions)

        //
        // annotate rhocall
        //
        vcf.map { meta, vcf, idx ->
                return [ vcf, idx ]
            }
            .set { ch_roh_vcfs}

        samples
            .combine(ch_roh_vcfs)
            .set { ch_roh_input }

        BCFTOOLS_ROH (ch_roh_input, gnomad_af, [], [], [], [])

    emit:
        vcf_ann                = VCFANNO.out.vcf
        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
