//
// Annotate with VEP
//

include { ENSEMBLVEP } from '../../modules/nf-core/modules/ensemblvep/main'

workflow ANNOTATE_VEP {
    take:
        vcf          // channel: [ val(meta), vcf ]
        vep_genome
        vep_species
        vep_cache_version
        vep_cache

    main:
        ch_reports  = Channel.empty()
        ch_vcf_ann  = Channel.empty()
        ch_versions = Channel.empty()

        ENSEMBLVEP(vcf, vep_genome, vep_species, vep_cache_version, vep_cache)

        ch_reports  = ch_reports.mix(ENSEMBLVEP.out.reports)
        ch_vcf_ann  = ch_vcf_ann.mix(ENSEMBLVEP.out.vcf_tbi)
        ch_versions = ch_versions.mix(ENSEMBLVEP.out.versions.first())

    emit:
        vcf_ann  = ch_vcf_ann      // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
        reports  = ch_reports      // path: *.html
        versions = ch_versions     // path: versions.yml}
