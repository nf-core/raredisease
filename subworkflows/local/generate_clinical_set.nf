//
// Generarte clinical set of variants
//

include { ENSEMBLVEP_FILTERVEP } from '../../modules/nf-core/ensemblvep/filtervep'
include { TABIX_BGZIP          } from '../../modules/nf-core/tabix/bgzip'
include { BCFTOOLS_FILTER      } from '../../modules/nf-core/bcftools/filter'

workflow GENERATE_CLINICAL_SET {
    take:
        ch_vcf      // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_hgnc_ids // channel: [mandatory] [ val(hgnc_ids) ]
        val_ismt    // value: if mitochondria, set to true

    main:
        ch_versions = Channel.empty()

        ENSEMBLVEP_FILTERVEP(
            ch_vcf,
            ch_hgnc_ids
        )
        .output
        .set { ch_filtervep_out }

        if (val_ismt) {
            BCFTOOLS_FILTER (ch_filtervep_out.map { meta, vcf -> return [meta, vcf, []]})
            ch_clinical = BCFTOOLS_FILTER.out.vcf
            ch_versions = ch_versions.mix( BCFTOOLS_FILTER.out.versions )
        } else {
            TABIX_BGZIP( ch_filtervep_out )
            ch_clinical = TABIX_BGZIP.out.output
            ch_versions = ch_versions.mix( TABIX_BGZIP.out.versions )
        }

        ch_versions = ch_versions.mix( ENSEMBLVEP_FILTERVEP.out.versions )

    emit:
        vcf      = ch_clinical // channel: [ val(meta), path(vcf) ]
        versions = ch_versions // channel: [ path(versions.yml) ]
}
