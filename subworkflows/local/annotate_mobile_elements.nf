//
// A subworkflow to annotate structural variants.
//

include { SVDB_QUERY as SVDB_QUERY_DB    } from '../../modules/nf-core/svdb/query/main'
include { PICARD_SORTVCF                 } from '../../modules/nf-core/picard/sortvcf/main'
include { ENSEMBLVEP as ENSEMBLVEP_ME    } from '../../modules/local/ensemblvep/main'
include { ENSEMBLVEP_FILTERVEP as FILTERVEP_ME  } from '../../modules/nf-core/ensemblvep/filtervep'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_FILTER } from '../../modules/nf-core/bcftools/view/main'
include { TABIX_BGZIPTABIX as BGZIP_TABIX_ME    } from '../../modules/nf-core/tabix/bgziptabix/main'

include { ANNOTATE_CSQ_PLI as ANNOTATE_CSQ_PLI_ME } from '../../subworkflows/local/annotate_consequence_pli.nf'

workflow ANNOTATE_MOBILE_ELEMENTS {

    take:
        ch_vcf                  // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_me_svdb_resources    // channel: [mandatory] [ path(csv) ]
        ch_genome_fasta         // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_dictionary    // channel: [mandatory] [ val(meta), path(dict) ]
        ch_vep_cache            // channel: [mandatory] [ path(cache) ]
        ch_variant_consequences // channel: [mandatory] [ path(consequences) ]
        ch_vep_filters          // channel: [mandatory] [ path(vep_filter) ]
        val_vep_genome          // string: [mandatory] GRCh37 or GRCh38
        val_vep_cache_version   // string: [mandatory] default: 107

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

        ENSEMBLVEP_ME(
            PICARD_SORTVCF.out.vcf,
            ch_genome_fasta,
            val_vep_genome,
            "homo_sapiens",
            val_vep_cache_version,
            ch_vep_cache,
            []
        )

        ENSEMBLVEP_ME.out.vcf_gz
            .map { meta, vcf ->
                [ meta, vcf, [] ]
            }
            .set { ch_bcftools_filter_input }

        BCFTOOLS_VIEW_FILTER ( ch_bcftools_filter_input, [], [], [] )

        BCFTOOLS_VIEW_FILTER.out.vcf
            .multiMap { meta, vcf ->
                clinical: [ meta + [ set: "clinical" ], vcf ]
                research: [ meta + [ set: "research" ], vcf ]   
            }
            .set { ch_clin_research_vcf }

        FILTERVEP_ME(
            ch_clin_research_vcf.clinical,
            ch_vep_filters
        )

        ch_clin_research_vcf.research
            .mix( FILTERVEP_ME.out.output )
            .set { ch_annotate_csq_pli_me_input }


        ANNOTATE_CSQ_PLI_ME (
            ch_annotate_csq_pli_me_input,
            ch_variant_consequences
        )

        ch_versions = ch_versions.mix(SVDB_QUERY_DB.out.versions)
        ch_versions = ch_versions.mix(PICARD_SORTVCF.out.versions)
        ch_versions = ch_versions.mix(ENSEMBLVEP_ME.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW_FILTER.out.versions)
        ch_versions = ch_versions.mix(FILTERVEP_ME.out.versions)
        ch_versions = ch_versions.mix(ANNOTATE_CSQ_PLI_ME.out.versions)

    emit:
        vcf      = ANNOTATE_CSQ_PLI_ME.out.vcf_ann  // channel: [ val(meta), path(vcf) ]
        tbi      = ANNOTATE_CSQ_PLI_ME.out.tbi_ann  // channel: [ val(meta), path(tbi) ]
        versions = ch_versions                      // channel: [ path(versions.yml) ]
}
