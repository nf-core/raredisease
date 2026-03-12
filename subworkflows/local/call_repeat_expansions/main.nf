//
// Run ExpansionHunter and Stranger
//

include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_EXP     } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_EXP   } from '../../../modules/nf-core/bcftools/reheader/main'
include { EXPANSIONHUNTER                              } from '../../../modules/nf-core/expansionhunter/main'
include { PICARD_RENAMESAMPLEINVCF as RENAMESAMPLE_EXP } from '../../../modules/nf-core/picard/renamesampleinvcf/main'
include { SAMTOOLS_SORT                                } from '../../../modules/nf-core/samtools/sort/main'
include { SVDB_MERGE as SVDB_MERGE_REPEATS             } from '../../../modules/nf-core/svdb/merge/main'
include { TABIX_TABIX as TABIX_EXP_RENAME              } from '../../../modules/nf-core/tabix/tabix/main'

workflow CALL_REPEAT_EXPANSIONS {
    take:
        ch_bam             // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_variant_catalog // channel: [mandatory] [ path(variant_catalog.json) ]
        ch_case_info       // channel: [mandatory] [ val(case_id) ]
        ch_genome_fasta    // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai      // channel: [mandatory] [ val(meta), path(fai) ]

    main:

        EXPANSIONHUNTER (
            ch_bam,
            ch_genome_fasta,
            ch_genome_fai,
            ch_variant_catalog
        )

        // Sort and index realigned bam
        SAMTOOLS_SORT(EXPANSIONHUNTER.out.bam, [[:],[]], 'bai')

        // Fix header and rename sample
        BCFTOOLS_REHEADER_EXP (
            EXPANSIONHUNTER.out.vcf.map{ meta, vcf -> [ meta, vcf, [], [] ]},
            ch_genome_fai
        )
        RENAMESAMPLE_EXP ( BCFTOOLS_REHEADER_EXP.out.vcf )
        TABIX_EXP_RENAME ( RENAMESAMPLE_EXP.out.vcf )

        // Split multi allelelic
        SPLIT_MULTIALLELICS_EXP (
            RENAMESAMPLE_EXP.out.vcf.join(TABIX_EXP_RENAME.out.index, failOnMismatch:true, failOnDuplicate:true),
            ch_genome_fasta
        )

        // Merge indiviual repeat expansions
        SPLIT_MULTIALLELICS_EXP.out.vcf
            .collect{_meta, vcf -> vcf}
            .toList()
            .collect()
            .set {ch_exp_vcfs}

        ch_case_info
            .combine(ch_exp_vcfs)
            .set {ch_svdb_merge_input}

        SVDB_MERGE_REPEATS ( ch_svdb_merge_input, [], true )

        ch_publish = SAMTOOLS_SORT.out.bam
            .mix(SAMTOOLS_SORT.out.bai)
            .mix(BCFTOOLS_REHEADER_EXP.out.vcf)
            .map { meta, value -> ['repeat_expansions/', [meta, value]] }

    emit:
        vcf        = SVDB_MERGE_REPEATS.out.vcf // channel: [ val(meta), path(vcf) ]
        ch_publish                              // channel: [ val(destination), val(value) ]
}
