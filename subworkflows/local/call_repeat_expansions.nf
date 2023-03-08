//
// Run ExpansionHunter and Stranger
//

include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_EXP     } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_EXP   } from '../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_VIEW as COMPRESS_STRANGER           } from '../../modules/nf-core/bcftools/view/main'
include { EXPANSIONHUNTER                              } from '../../modules/nf-core/expansionhunter/main'
include { PICARD_RENAMESAMPLEINVCF as RENAMESAMPLE_EXP } from '../../modules/nf-core/picard/renamesampleinvcf/main'
include { STRANGER                                     } from '../../modules/nf-core/stranger/main'
include { SVDB_MERGE as SVDB_MERGE_REPEATS             } from '../../modules/nf-core/svdb/merge/main'
include { TABIX_BGZIPTABIX as BGZIPTABIX_EXP           } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX as INDEX_STRANGER                } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_EXP_RENAME              } from '../../modules/nf-core/tabix/tabix/main'

workflow CALL_REPEAT_EXPANSIONS {
    take:
        bam 			// channel: [ val(meta), path(bam), path(bai) ]
        variant_catalog	// channel: /path/to/variant_catalog.json
        case_info       // channel: [ case_id ]
        fasta 			// channel: /path/to/genome.fasta
        fai 			// channel: /path/to/genome.fasta.fai

    main:
        ch_versions = Channel.empty()

        EXPANSIONHUNTER (
            bam,
            fasta,
            variant_catalog
        )

        // Fix header and rename sample
        BCFTOOLS_REHEADER_EXP (
            EXPANSIONHUNTER.out.vcf.map{ meta, vcf -> [ meta, vcf, [] ]},
            fai
        )
        RENAMESAMPLE_EXP ( BCFTOOLS_REHEADER_EXP.out.vcf )
        TABIX_EXP_RENAME ( RENAMESAMPLE_EXP.out.vcf )

        // Split multi allelelic
        SPLIT_MULTIALLELICS_EXP (
            RENAMESAMPLE_EXP.out.vcf.join(TABIX_EXP_RENAME.out.tbi),
            fasta
        )

        // Merge indiviual repeat expansions
        SPLIT_MULTIALLELICS_EXP.out.vcf
            .collect{it[1]}
            .toList()
            .set {ch_exp_vcfs}
        case_info
            .combine(ch_exp_vcfs)
            .set {ch_svdb_merge_input}
        SVDB_MERGE_REPEATS ( ch_svdb_merge_input, [] )

        // Annotate, compress and index
        STRANGER ( SVDB_MERGE_REPEATS.out.vcf, variant_catalog )
        COMPRESS_STRANGER (
            STRANGER.out.vcf.map{ meta, vcf -> [meta, vcf, [] ]},
             [], [], []
        )
        INDEX_STRANGER ( COMPRESS_STRANGER.out.vcf )

        ch_versions = ch_versions.mix(EXPANSIONHUNTER.out.versions.first())
        ch_versions = ch_versions.mix(BCFTOOLS_REHEADER_EXP.out.versions.first())
        ch_versions = ch_versions.mix(RENAMESAMPLE_EXP.out.versions.first()    )
        ch_versions = ch_versions.mix(TABIX_EXP_RENAME.out.versions.first())
        ch_versions = ch_versions.mix(SPLIT_MULTIALLELICS_EXP.out.versions.first())
        ch_versions = ch_versions.mix(SVDB_MERGE_REPEATS.out.versions.first())
        ch_versions = ch_versions.mix(STRANGER.out.versions.first())
        ch_versions = ch_versions.mix(COMPRESS_STRANGER.out.versions.first())
        ch_versions = ch_versions.mix(INDEX_STRANGER.out.versions.first())

    emit:
        vcf         = COMPRESS_STRANGER.out.vcf.join(INDEX_STRANGER.out.tbi) // channel: [ val(meta), path(*.vcf.gz), path(*.vcf.gz.tbi) ]
        versions    = ch_versions                                            // channel: [ versions.yml ]
}
