//
// A subworkflow to annotate structural variants.
//

include { SVDB_QUERY as SVDB_QUERY_FILTER } from '../../modules/nf-core/modules/svdb/query/main'
include { SVDB_QUERY as SVDB_QUERY_NOFILTER } from '../../modules/nf-core/modules/svdb/query/main'


// include { PICARD_SORTVCF } from '../../modules/nf-core/modules/picard/sortvcf/main'

workflow ANNOTATE_STRUCTURAL_VARIANTS {

    take:
        vcf         // channel: [ val(meta), path(vcf) ]
        sv_dbs      // file: dbs.csv

    main:
        ch_versions = Channel.empty()

        Channel.fromPath(sv_dbs)
            .splitCsv ( header:true )
            .branch { row ->
                freq_filter: row.use_in_freq_filter == "1"
                    return [row.filename,
                            row.in_freq_info_key,
                            row.in_allele_count_info_key,
                            row.out_freq_info_key,
                            row.out_allele_count_info_key]
                no_freq_filter: row.use_in_freq_filter == "0"
                    return [row.filename,
                            row.in_freq_info_key,
                            row.in_allele_count_info_key,
                            row.out_freq_info_key,
                            row.out_allele_count_info_key]
            }
            .set { ch_svdb_dbs }

        ch_input_filter   = vcf.combine(ch_svdb_dbs.freq_filter).view()
        ch_input_nofilter = vcf.combine(ch_svdb_dbs.no_freq_filter).view()

        SVDB_QUERY_FILTER(ch_input_filter)
        SVDB_QUERY_NOFILTER(ch_input_nofilter)

    emit:
        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
