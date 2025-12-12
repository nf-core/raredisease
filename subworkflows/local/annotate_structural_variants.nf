//
// A subworkflow to annotate structural variants.
//

include { SVDB_QUERY as SVDB_QUERY_DB     } from '../../modules/nf-core/svdb/query/main'
include { SVDB_QUERY as SVDB_QUERY_BEDPE  } from '../../modules/nf-core/svdb/query/main'
include { PICARD_SORTVCF                  } from '../../modules/nf-core/picard/sortvcf/main'
include { BCFTOOLS_VIEW                   } from '../../modules/nf-core/bcftools/view/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_SV } from '../../modules/nf-core/ensemblvep/vep/main'
include { TABIX_TABIX as TABIX_VEP        } from '../../modules/nf-core/tabix/tabix/main'

workflow ANNOTATE_STRUCTURAL_VARIANTS {

    take:
        ch_vcf                // channel: [mandatory] [ val(meta), path(vcf) ]
        svdb_query_bedpedbs   // String: [optional] params.svdb_query_bedpedbs
        svdb_query_dbs        // String: [optional] params.svdb_query_dbs
        val_vep_genome        // string: [mandatory] GRCh37 or GRCh38
        val_vep_cache_version // string: [mandatory] default: 107
        ch_svdb_bedpedbs      // channel: [optional] 
        ch_svdb_dbs           // channel: [optional] 
        ch_vep_cache          // channel: [mandatory] [ path(cache) ]
        ch_genome_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_dictionary  // channel: [mandatory] [ val(meta), path(dict) ]
        ch_vep_extra_files    // channel: [mandatory] [ path(files) ]

    main:
        ch_versions      = channel.empty()

        if (svdb_query_dbs) {
            ch_svdb_dbs
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

            ch_vcf      = SVDB_QUERY_DB.out.vcf
            ch_versions = ch_versions.mix(SVDB_QUERY_DB.out.versions)
        }

        if (svdb_query_bedpedbs) {
            ch_svdb_bedpedbs
                .multiMap { file, in_freq_info_key, in_allele_count_info_key, out_freq_info_key, out_allele_count_info_key ->
                    bedpedbs: file
                    in_frqs:  in_freq_info_key
                    in_occs:  in_allele_count_info_key
                    out_frqs: out_freq_info_key
                    out_occs: out_allele_count_info_key
                }
                .set { ch_svdb_bedpedbs }

            SVDB_QUERY_BEDPE (
                ch_vcf,
                ch_svdb_bedpedbs.in_occs.toList(),
                ch_svdb_bedpedbs.in_frqs.toList(),
                ch_svdb_bedpedbs.out_occs.toList(),
                ch_svdb_bedpedbs.out_frqs.toList(),
                [],
                ch_svdb_bedpedbs.bedpedbs.toList()
            )

            ch_vcf      = SVDB_QUERY_BEDPE.out.vcf
            ch_versions = ch_versions.mix(SVDB_QUERY_BEDPE.out.versions)
        }

        PICARD_SORTVCF(ch_vcf, ch_genome_fasta, ch_genome_dictionary)

        PICARD_SORTVCF.out.vcf
            .map { meta, vcf -> return [meta,vcf,[]] }
            .set { ch_sortvcf }

        BCFTOOLS_VIEW(ch_sortvcf, [], [], [])
            .vcf
            .map { meta, vcf -> return [meta, vcf, []]}
            .set { ch_vep_in }

        ENSEMBLVEP_SV(
            ch_vep_in,
            val_vep_genome,
            "homo_sapiens",
            val_vep_cache_version,
            ch_vep_cache,
            ch_genome_fasta,
            ch_vep_extra_files
        )

        TABIX_VEP (ENSEMBLVEP_SV.out.vcf)

        ch_versions = ch_versions.mix(PICARD_SORTVCF.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)
        ch_versions = ch_versions.mix(ENSEMBLVEP_SV.out.versions)
        ch_versions = ch_versions.mix(TABIX_VEP.out.versions)

    emit:
        vcf_ann  = ENSEMBLVEP_SV.out.vcf // channel: [ val(meta), path(vcf) ]
        tbi      = TABIX_VEP.out.tbi     // channel: [ val(meta), path(tbi) ]
        versions = ch_versions           // channel: [ path(versions.yml) ]
}
