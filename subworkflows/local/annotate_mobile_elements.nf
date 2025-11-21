//
// A subworkflow to annotate structural variants.
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_FILTER } from '../../modules/nf-core/bcftools/view/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_ME       } from '../../modules/nf-core/ensemblvep/vep/main'
include { PICARD_SORTVCF                        } from '../../modules/nf-core/picard/sortvcf/main'
include { SVDB_QUERY as SVDB_QUERY_DB           } from '../../modules/nf-core/svdb/query/main'


workflow ANNOTATE_MOBILE_ELEMENTS {

    take:
        ch_vcf                  // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_me_svdb_resources    // channel: [mandatory] [ path(csv) ]
        ch_genome_fasta         // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_dictionary    // channel: [mandatory] [ val(meta), path(dict) ]
        ch_vep_cache            // channel: [mandatory] [ path(cache) ]
        val_vep_genome          // string: [mandatory] GRCh37 or GRCh38
        val_vep_cache_version   // string: [mandatory] default: 107
        ch_vep_extra_files      // channel: [mandatory] [ path(files) ]

    main:
        ch_versions = Channel.empty()
        ch_svdb_dbs = Channel.empty()

        ch_me_svdb_resources
            .splitCsv ( header:true )
            .multiMap { row ->
                vcf_dbs:  row.filename
                in_frqs:  row.in_freq_info_key
                in_occs:  row.in_allele_count_info_key
                out_frqs: row.out_freq_info_key
                out_occs: row.out_allele_count_info_key
            }
            .set { ch_svdb_dbs }

        SVDB_QUERY_DB (
            ch_vcf,
            ch_svdb_dbs.in_occs.toList(),
            ch_svdb_dbs.in_frqs.toList(),
            ch_svdb_dbs.out_occs.toList(),
            ch_svdb_dbs.out_frqs.toList(),
            ch_svdb_dbs.vcf_dbs.toList(),
            []
        )

        PICARD_SORTVCF(
            SVDB_QUERY_DB.out.vcf,
            ch_genome_fasta,
            ch_genome_dictionary
        )
        .vcf
        .map { meta, vcf -> return [meta, vcf, []] }
        .set { ch_vep_in }

        ENSEMBLVEP_ME(
            ch_vep_in,
            val_vep_genome,
            "homo_sapiens",
            val_vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta,
            ch_vep_extra_files
        )

        ENSEMBLVEP_ME.out.vcf
            .map { meta, vcf ->
                [ meta, vcf, [] ]
            }
            .set { ch_bcftools_filter_input }

        BCFTOOLS_VIEW_FILTER( ch_bcftools_filter_input, [], [], [] )

        ch_versions = ch_versions.mix( SVDB_QUERY_DB.out.versions )
        ch_versions = ch_versions.mix( PICARD_SORTVCF.out.versions )
        ch_versions = ch_versions.mix( ENSEMBLVEP_ME.out.versions )
        ch_versions = ch_versions.mix( BCFTOOLS_VIEW_FILTER.out.versions )

    emit:
        vcf_ann  = BCFTOOLS_VIEW_FILTER.out.vcf     // channel: [ val(meta), path(vcf) ]
        versions = ch_versions                      // channel: [ path(versions.yml) ]
}
