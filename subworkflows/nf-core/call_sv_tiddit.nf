//
// A structural variant caller workflow for tiddit
//

include { TIDDIT_SV } from '../../modules/nf-core/tiddit/sv/main'

include { SVDB_MERGE as SVDB_MERGE_TIDDIT } from '../../modules/nf-core/svdb/merge/main'

workflow CALL_SV_TIDDIT {
    take:
    bam            // channel: [ val(meta), path(bam) ]
    fasta          // path(fasta)
    index          // [ val(meta), path(bwa_index)]
    case_info      // channel: [ case_id ]

    main:
        index_for_tiddit = index.map { meta, ind -> ind }
        TIDDIT_SV ( bam, fasta, index_for_tiddit )
        ch_versions = TIDDIT_SV.out.versions

        TIDDIT_SV.out
            .vcf
            .collect{it[1]}
            .toList()
            .set { vcf_file_list }

        case_info
            .combine(vcf_file_list)
            .set { merge_input_vcfs }

        SVDB_MERGE_TIDDIT ( merge_input_vcfs, [] )
        ch_versions = ch_versions.mix(SVDB_MERGE_TIDDIT.out.versions)

    emit:
        vcf          = SVDB_MERGE_TIDDIT.out.vcf
        versions     = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
