//
// A variant caller workflow for deepvariant
//

include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_GL } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES_GL   } from '../../../modules/nf-core/bcftools/norm/main'
include { DEEPVARIANT                             } from '../../../modules/nf-core/deepvariant/main'
include { GLNEXUS                                 } from '../../../modules/nf-core/glnexus/main'
include { TABIX_TABIX as TABIX_GL                 } from '../../../modules/nf-core/tabix/tabix/main'

workflow CALL_SNV_DEEPVARIANT {
    take:
        ch_bam_bai      // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai   // channel: [mandatory] [ val(meta), path(fai) ]
        ch_case_info    // channel: [mandatory] [ val(case_info) ]

    main:
        ch_versions = Channel.empty()

        ch_bam_bai.map { meta, bam, bai ->
                        return [meta, bam, bai, []]
            }
            .set { ch_deepvar_in }

        DEEPVARIANT ( ch_deepvar_in, ch_genome_fasta, ch_genome_fai, [[],[]] )
        DEEPVARIANT.out.gvcf
            .collect{it[1]}
            .toList()
            .collect()
            .set { ch_file_list }

        ch_case_info
            .combine(ch_file_list)
            .set { ch_gvcfs }

        GLNEXUS ( ch_gvcfs )

        ch_split_multi_in = GLNEXUS.out.bcf
                            .map{ meta, bcf ->
                                    return [meta, bcf, []] }
        SPLIT_MULTIALLELICS_GL (ch_split_multi_in, ch_genome_fasta)

        ch_remove_dup_in = SPLIT_MULTIALLELICS_GL.out.vcf
                            .map{ meta, vcf ->
                                    return [meta, vcf, []] }
        REMOVE_DUPLICATES_GL (ch_remove_dup_in, ch_genome_fasta)

        TABIX_GL (REMOVE_DUPLICATES_GL.out.vcf)

        ch_versions = ch_versions.mix(DEEPVARIANT.out.versions.first())
        ch_versions = ch_versions.mix(GLNEXUS.out.versions)
        ch_versions = ch_versions.mix(SPLIT_MULTIALLELICS_GL.out.versions)
        ch_versions = ch_versions.mix(REMOVE_DUPLICATES_GL.out.versions)
        ch_versions = ch_versions.mix(TABIX_GL.out.versions)

    emit:
        vcf      = REMOVE_DUPLICATES_GL.out.vcf // channel: [ val(meta), path(vcf) ]
        tabix    = TABIX_GL.out.tbi             // channel: [ val(meta), path(tbi) ]
        versions = ch_versions                  // channel: [ path(versions.yml) ]
}
