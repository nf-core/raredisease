//
// A structural variant caller workflow for tiddit
//

include { TIDDIT_SV                       } from '../../../modules/nf-core/tiddit/sv/main'
include { SVDB_MERGE as SVDB_MERGE_TIDDIT } from '../../../modules/nf-core/svdb/merge/main'

workflow CALL_SV_TIDDIT {
    take:
        ch_bam_bai      // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_bwa_index    // channel: [mandatory] [ val(meta), path(index)]
        ch_case_info    // channel: [mandatory] [ val(case_info) ]

    main:
        TIDDIT_SV ( ch_bam_bai, ch_genome_fasta, ch_bwa_index )

        TIDDIT_SV.out
            .vcf
            .collect{it[1]}
            .toList()
            .set { vcf_file_list }

        ch_case_info
            .combine(vcf_file_list)
            .set { merge_input_vcfs }

        SVDB_MERGE_TIDDIT ( merge_input_vcfs, [] )

        ch_versions = TIDDIT_SV.out.versions.first()
        ch_versions = ch_versions.mix(SVDB_MERGE_TIDDIT.out.versions)

    emit:
        vcf      = SVDB_MERGE_TIDDIT.out.vcf // channel: [ val(meta), path(vcf) ]
        versions = ch_versions               // channel: [ path(versions.yml) ]
}
