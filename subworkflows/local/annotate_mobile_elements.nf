//
// A subworkflow to annotate structural variants.
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_FILTER } from '../../modules/nf-core/bcftools/view/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_ME       } from '../../modules/nf-core/ensemblvep/vep/main'
include { PICARD_SORTVCF                        } from '../../modules/nf-core/picard/sortvcf/main'
include { SVDB_QUERY as SVDB_QUERY_DB           } from '../../modules/nf-core/svdb/query/main'


workflow ANNOTATE_MOBILE_ELEMENTS {

    take:
        ch_genome_dictionary    // channel: [mandatory] [ val(meta), path(dict) ]
        ch_genome_fasta         // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_me_svdb_resources    // channel: [mandatory] [ path(csv) ]
        ch_vcf                  // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_vep_cache            // channel: [mandatory] [ path(cache) ]
        val_genome              // string: [mandatory] GRCh37 or GRCh38
        val_vep_cache_version   // string: [mandatory] default: 107
        ch_vep_extra_files      // channel: [mandatory] [ path(files) ]

    main:
        ch_svdb_dbs = channel.empty()

        ch_me_svdb_resources
            .multiMap { file, in_freq_info_key, in_allele_count_info_key, out_freq_info_key, out_allele_count_info_key ->
                vcf_dbs:  file
                in_frqs:  in_freq_info_key
                in_occs:  in_allele_count_info_key
                out_frqs: out_freq_info_key
                out_occs: out_allele_count_info_key
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
            val_genome,
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

    emit:
        vcf_ann  = BCFTOOLS_VIEW_FILTER.out.vcf     // channel: [ val(meta), path(vcf) ]
}
