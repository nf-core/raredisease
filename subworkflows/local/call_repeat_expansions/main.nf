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
        ch_samplegender    // channel: [mandatory] [ val(meta), path(tsv) ]
        ch_variant_catalog // channel: [mandatory] [ path(variant_catalog.json) ]
        ch_case_info       // channel: [mandatory] [ val(case_id) ]
        ch_genome_fasta    // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai      // channel: [mandatory] [ val(meta), path(fai) ]

    main:
        ch_samplegender_parsed = ch_samplegender.map { meta, tsv ->
            def data_line = tsv.readLines()
                .find { line -> line.trim() && !line.startsWith('#') }

            def gender
            if (data_line) {
                gender = data_line.split('\t')[1]
            } else if (workflow.stubRun) {
                gender = meta.sex?.toString() == '2' ? 'female' : 'male'
            } else {
                throw new IllegalStateException(
                    "No SampleGender result found for sample ${meta.id} in ${tsv}"
                )
            }

            [meta.id, gender]
        }

       ch_bam_with_gender = ch_bam
           .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
           .join(ch_samplegender_parsed, by: 0, failOnMismatch: true, failOnDuplicate: true)
           .map { id, meta, bam, bai, gender ->
               [meta + [ngsbits_sex: gender], bam, bai]
       }

         EXPANSIONHUNTER (
            ch_bam_with_gender,
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
        ch_exp_vcfs = SPLIT_MULTIALLELICS_EXP.out.vcf
            .collect{_meta, vcf -> vcf}
            .toList()
            .collect()

        ch_svdb_merge_input = ch_case_info
            .combine(ch_exp_vcfs)

        SVDB_MERGE_REPEATS ( ch_svdb_merge_input, [], true )

    emit:
        expansionhunter_bai = SAMTOOLS_SORT.out.bai           // channel: [ val(meta), path(bai) ]
        expansionhunter_bam = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), path(bam) ]
        expansionhunter_vcf = BCFTOOLS_REHEADER_EXP.out.vcf  // channel: [ val(meta), path(vcf) ]
        vcf                 = SVDB_MERGE_REPEATS.out.vcf      // channel: [ val(meta), path(vcf) ]
}
