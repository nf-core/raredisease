//
// A structural variant caller workflow for tiddit
//

include { TIDDIT_SV                             } from '../../../modules/nf-core/tiddit/sv/main'
include { SVDB_MERGE as SVDB_MERGE_TIDDIT       } from '../../../modules/nf-core/svdb/merge/main'
include { TABIX_BGZIPTABIX as INDEX_TIDDIT      } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_TIDDIT } from '../../../modules/nf-core/bcftools/view/main.nf'

workflow CALL_SV_TIDDIT {
    take:
        ch_bam_bai      // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_fai   // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_bwa_index    // channel: [mandatory] [ val(meta), path(index)]
        ch_case_info    // channel: [mandatory] [ val(case_info) ]

    main:
        ch_genome_fa_fai = ch_genome_fasta.join(ch_genome_fai).collect()
        TIDDIT_SV ( ch_bam_bai, ch_genome_fa_fai, ch_bwa_index )

        INDEX_TIDDIT (TIDDIT_SV.out.vcf)
        BCFTOOLS_VIEW_TIDDIT (INDEX_TIDDIT.out.gz_index, [], [], []).vcf
            .collect{ _meta, vcf -> vcf}
            .toList()
            .set { vcf_file_list }

        ch_case_info
            .combine(vcf_file_list)
            .set { merge_input_vcfs }

        SVDB_MERGE_TIDDIT ( merge_input_vcfs, [], true )

    emit:
        vcf      = SVDB_MERGE_TIDDIT.out.vcf // channel: [ val(meta), path(vcf) ]
}
