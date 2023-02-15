//
// A subworkflow to annotate structural variants.
//

include { SVDB_QUERY                    } from '../../modules/nf-core/svdb/query/main'
include { PICARD_SORTVCF                } from '../../modules/nf-core/picard/sortvcf/main'
include { BCFTOOLS_VIEW                 } from '../../modules/nf-core/bcftools/view/main'
include { ENSEMBLVEP as ENSEMBLVEP_SV   } from '../../modules/local/ensemblvep/main'

workflow ANNOTATE_STRUCTURAL_VARIANTS {

    take:
        vcf               // channel: [ val(meta), path(vcf) ]
        sv_dbs            // file: dbs.csv
        vep_genome
        vep_cache_version
        vep_cache
        fasta             // file: genome.fasta
        seq_dict          // file: genome.dict

    main:
        ch_versions = Channel.empty()

        Channel.fromPath(sv_dbs)
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
            vcf,
            ch_svdb_dbs.in_occs.toList(),
            ch_svdb_dbs.in_frqs.toList(),
            ch_svdb_dbs.out_occs.toList(),
            ch_svdb_dbs.out_frqs.toList(),
            ch_svdb_dbs.vcf_dbs.toList()
        )

        PICARD_SORTVCF(SVDB_QUERY.out.vcf, fasta, seq_dict)

        PICARD_SORTVCF.out.vcf.map { meta, vcf -> return [meta,vcf,[]] }.set { ch_sortvcf }

        BCFTOOLS_VIEW(ch_sortvcf, [], [], [])

        ENSEMBLVEP_SV(
            BCFTOOLS_VIEW.out.vcf,
            vep_genome,
            "homo_sapiens",
            vep_cache_version,
            vep_cache,
            fasta,
            []
        )

        ch_versions = ch_versions.mix(SVDB_QUERY.out.versions)
        ch_versions = ch_versions.mix(PICARD_SORTVCF.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)
        ch_versions = ch_versions.mix(ENSEMBLVEP_SV.out.versions)

    emit:
        vcf_ann       = ENSEMBLVEP_SV.out.vcf
        versions      = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
