//
// A variant caller workflow for deepvariant
//

include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_GL } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES_GL   } from '../../modules/nf-core/bcftools/norm/main'
include { DEEPVARIANT                             } from '../../modules/nf-core/deepvariant/main'
include { GLNEXUS                                 } from '../../modules/nf-core/glnexus/main'
include { TABIX_TABIX as TABIX_GL                 } from '../../modules/nf-core/tabix/tabix/main'

workflow CALL_SNV_DEEPVARIANT {
    take:
        bam          // channel: [ val(meta), path(bam), path(bai) ]
        fasta        // path(fasta)
        fai          // path(fai)
        case_info    // channel: [ case_id ]

    main:
        ch_versions = Channel.empty()
        bam.map { meta, bam, bai ->
                        return [meta, bam, bai, []]
            }
            .set { ch_bam }

        DEEPVARIANT ( ch_bam, fasta, fai )
        DEEPVARIANT.out
            .gvcf
            .collect{it[1]}
            .toList()
            .set { file_list }
        ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)

        case_info
            .combine(file_list)
            .set { ch_gvcfs }

        GLNEXUS ( ch_gvcfs )
        ch_versions = ch_versions.mix(GLNEXUS.out.versions)

        ch_split_multi_in = GLNEXUS.out.bcf
                            .map{meta, bcf ->
                                    return [meta, bcf, []]}
        SPLIT_MULTIALLELICS_GL (ch_split_multi_in, fasta)
        ch_versions = ch_versions.mix(SPLIT_MULTIALLELICS_GL.out.versions)

        ch_remove_dup_in = SPLIT_MULTIALLELICS_GL.out.vcf
                            .map{meta, vcf ->
                                    return [meta, vcf, []]}
        REMOVE_DUPLICATES_GL (ch_remove_dup_in, fasta)
        ch_versions = ch_versions.mix(REMOVE_DUPLICATES_GL.out.versions)

        TABIX_GL (REMOVE_DUPLICATES_GL.out.vcf)
        ch_versions = ch_versions.mix(TABIX_GL.out.versions)

    emit:
        vcf         = REMOVE_DUPLICATES_GL.out.vcf
        tabix       = TABIX_GL.out.tbi
        versions    = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
