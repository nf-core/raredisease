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
        ch_sv_dbs             // channel: [mandatory] [ val(csv) ]
        ch_sv_bedpedbs        // channel: [mandatory] [ val(csv) ]
        val_vep_genome        // string: [mandatory] GRCh37 or GRCh38
        val_vep_cache_version // string: [mandatory] default: 107
        ch_vep_cache          // channel: [mandatory] [ path(cache) ]
        ch_genome_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_dictionary  // channel: [mandatory] [ val(meta), path(dict) ]
        ch_vep_extra_files    // channel: [mandatory] [ path(files) ]

    main:
        ch_versions      = Channel.empty()
        ch_svdb_dbs      = Channel.empty()
        ch_svdb_bedpedbs = Channel.empty()

        ch_sv_dbs
            .splitCsv ( header:true )
            .multiMap { row ->
                vcf_dbs:  row.filename
                in_frqs:  row.in_freq_info_key
                in_occs:  row.in_allele_count_info_key
                out_frqs: row.out_freq_info_key
                out_occs: row.out_allele_count_info_key
            }
            .set { ch_svdb_dbs }

        ch_sv_bedpedbs
            .splitCsv ( header:true )
            .multiMap { row ->
                bedpedbs: row.filename
                in_frqs:  row.in_freq_info_key
                in_occs:  row.in_allele_count_info_key
                out_frqs: row.out_freq_info_key
                out_occs: row.out_allele_count_info_key
            }
            .set { ch_svdb_bedpedbs }

        SVDB_QUERY_DB (
            ch_vcf,
            ch_svdb_dbs.in_occs.toList(),
            ch_svdb_dbs.in_frqs.toList(),
            ch_svdb_dbs.out_occs.toList(),
            ch_svdb_dbs.out_frqs.toList(),
            ch_svdb_dbs.vcf_dbs.toList(),
            []
        )

        ch_vcf
            .join(SVDB_QUERY_DB.out.vcf, remainder: true)
            .branch { it  ->
                original_call: it[2].equals(null)
                    return [it[0], it[1]]
                annotated_with_db: !(it[2].equals(null))
                    return [it[0], it[2]]
            }
            .set { ch_for_mix_querydb }

        ch_querydb_out = ch_for_mix_querydb.original_call.mix(ch_for_mix_querydb.annotated_with_db)

        SVDB_QUERY_BEDPE (
            ch_querydb_out,
            ch_svdb_bedpedbs.in_occs.toList(),
            ch_svdb_bedpedbs.in_frqs.toList(),
            ch_svdb_bedpedbs.out_occs.toList(),
            ch_svdb_bedpedbs.out_frqs.toList(),
            [],
            ch_svdb_bedpedbs.bedpedbs.toList()
        )

        ch_querydb_out
            .join(SVDB_QUERY_BEDPE.out.vcf, remainder: true)
            .branch { it  ->
                querydb_out: it[2].equals(null)
                    return [it[0], it[1]]
                annotated_with_bedped: !(it[2].equals(null))
                    return [it[0], it[2]]
            }
            .set { ch_for_mix_querybedpedb }

        ch_querypedbed_out = ch_for_mix_querybedpedb.querydb_out.mix(ch_for_mix_querybedpedb.annotated_with_bedped)

        PICARD_SORTVCF(ch_querypedbed_out, ch_genome_fasta, ch_genome_dictionary)

        PICARD_SORTVCF.out.vcf.map { meta, vcf -> return [meta,vcf,[]] }.set { ch_sortvcf }

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

        ch_versions = ch_versions.mix(SVDB_QUERY_DB.out.versions)
        ch_versions = ch_versions.mix(SVDB_QUERY_BEDPE.out.versions)
        ch_versions = ch_versions.mix(PICARD_SORTVCF.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)
        ch_versions = ch_versions.mix(ENSEMBLVEP_SV.out.versions)
        ch_versions = ch_versions.mix(TABIX_VEP.out.versions)

    emit:
        vcf_ann  = ENSEMBLVEP_SV.out.vcf // channel: [ val(meta), path(vcf) ]
        tbi      = TABIX_VEP.out.tbi     // channel: [ val(meta), path(tbi) ]
        versions = ch_versions           // channel: [ path(versions.yml) ]
}
