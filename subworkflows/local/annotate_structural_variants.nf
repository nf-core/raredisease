//
// A subworkflow to annotate structural variants.
//

include { SVDB_QUERY                  } from '../../modules/nf-core/svdb/query/main'
include { PICARD_SORTVCF              } from '../../modules/nf-core/picard/sortvcf/main'
include { BCFTOOLS_VIEW               } from '../../modules/nf-core/bcftools/view/main'
include { ENSEMBLVEP as ENSEMBLVEP_SV } from '../../modules/local/ensemblvep/main'
include { TABIX_TABIX as TABIX_VEP    } from '../../modules/nf-core/tabix/tabix/main'

workflow ANNOTATE_STRUCTURAL_VARIANTS {

    take:
        ch_vcf                // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_sv_dbs             // channel: [mandatory] [ val(csv) ]
        val_vep_genome        // string: [mandatory] GRCh37 or GRCh38
        val_vep_cache_version // string: [mandatory] default: 107
        ch_vep_cache          // channel: [mandatory] [ path(cache) ]
        ch_genome_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_dictionary  // channel: [mandatory] [ val(meta), path(dict) ]

    main:
        ch_versions = Channel.empty()

        Channel.fromPath(ch_sv_dbs)
            .splitCsv ( header:true )
            .multiMap { row ->
                vcf_dbs:  row.filename
                in_frqs:  row.in_freq_info_key
                in_occs:  row.in_allele_count_info_key
                out_frqs: row.out_freq_info_key
                out_occs: row.out_allele_count_info_key
            }
            .set { ch_svdb_dbs }

        SVDB_QUERY(
            ch_vcf,
            ch_svdb_dbs.in_occs.toList(),
            ch_svdb_dbs.in_frqs.toList(),
            ch_svdb_dbs.out_occs.toList(),
            ch_svdb_dbs.out_frqs.toList(),
            ch_svdb_dbs.vcf_dbs.toList()
        )

        PICARD_SORTVCF(SVDB_QUERY.out.vcf, ch_genome_fasta, ch_genome_dictionary)

        PICARD_SORTVCF.out.vcf.map { meta, vcf -> return [meta,vcf,[]] }.set { ch_sortvcf }

        BCFTOOLS_VIEW(ch_sortvcf, [], [], [])

        ENSEMBLVEP_SV(
            BCFTOOLS_VIEW.out.vcf,
            ch_genome_fasta,
            val_vep_genome,
            "homo_sapiens",
            val_vep_cache_version,
            ch_vep_cache,
            []
        )

        TABIX_VEP (ENSEMBLVEP_SV.out.vcf_gz)

        ch_versions = ch_versions.mix(SVDB_QUERY.out.versions)
        ch_versions = ch_versions.mix(PICARD_SORTVCF.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)
        ch_versions = ch_versions.mix(ENSEMBLVEP_SV.out.versions)
        ch_versions = ch_versions.mix(TABIX_VEP.out.versions)

    emit:
        vcf_ann  = ENSEMBLVEP_SV.out.vcf_gz // channel: [ val(meta), path(vcf) ]
        tbi      = TABIX_VEP.out.tbi        // channel: [ val(meta), path(tbi) ]
        versions = ch_versions              // channel: [ path(versions.yml) ]
}
