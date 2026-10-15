//
// A subworkflow to call mobile elements in the genome
//

include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_ME     } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_ME } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_ME         } from '../../../modules/nf-core/bcftools/sort/main'
include { RETROSEQ_CALL as RETROSEQ_CALL            } from '../../../modules/local/retroseq/call/main'
include { RETROSEQ_DISCOVER as RETROSEQ_DISCOVER    } from '../../../modules/local/retroseq/discover/main'
include { SAMTOOLS_VIEW as ME_SPLIT_ALIGNMENT       } from '../../../modules/nf-core/samtools/view/main'
include { SVDB_MERGE as SVDB_MERGE_ME               } from '../../../modules/nf-core/svdb/merge/main'
include { TABIX_TABIX as TABIX_ME                   } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_ME_SPLIT             } from '../../../modules/nf-core/tabix/tabix/main'

workflow CALL_MOBILE_ELEMENTS {

    take:
        ch_case_info        // channel: [mandatory] [ val(case_info) ]
        ch_genome_bam_bai   // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_fai       // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_fasta     // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_me_references    // channel: [mandatory] [path(tsv)]

    main:
        // Building chromosome channels based on fasta index
        ch_chr = ch_genome_fai
            .splitCsv( sep: "\t", elem: 1, limit: 24 )
            .map { _meta, fai -> [ fai.first() ] }
            .collect()
            .map { chr -> [ chr, chr.size() ] }
            .transpose()

        // Building one bam channel per chromosome and adding interval and the number of intervals
        ch_genome_bam_bai_interval = ch_genome_bam_bai
            .combine( ch_chr )
            .map { meta, bam, bai, chr, nr_of_chrs ->
                [ meta + [interval:chr, nr_of_intervals: nr_of_chrs], bam, bai ]
            }

        // Split bam file on chromosome and index
        ME_SPLIT_ALIGNMENT ( ch_genome_bam_bai_interval, [[:], [], []], [], 'bai' )

        ch_retroseq_input = ME_SPLIT_ALIGNMENT.out.bam
            .join( ME_SPLIT_ALIGNMENT.out.bai, failOnMismatch: true, failOnDuplicate: true )

        ch_me_reference_split = ch_me_references
            .multiMap { type, path ->
                type: type
                path: path
            }

        RETROSEQ_DISCOVER (
            ch_retroseq_input,
            ch_me_reference_split.path.collect(),
            ch_me_reference_split.type.collect()
        )

        ch_retroseq_call_input = RETROSEQ_DISCOVER.out.tab
            .join(ch_retroseq_input, failOnMismatch: true)

        RETROSEQ_CALL (
            ch_retroseq_call_input,
            ch_genome_fasta,
            ch_genome_fai
        )

        // Fix the vcf by adding header, sorting and indexing
        BCFTOOLS_REHEADER_ME (
            RETROSEQ_CALL.out.vcf.map{ meta, vcf -> [ meta, vcf, [], [] ] },
            ch_genome_fai
        )
        BCFTOOLS_SORT_ME ( BCFTOOLS_REHEADER_ME.out.vcf )
        TABIX_ME_SPLIT ( BCFTOOLS_SORT_ME.out.vcf )

        // Preparing channels for input to bcftools concat
        // resulting channel [ meta, [ vcf_1, vcf_2, ... ], [ tbi_1, tbi_2, ... ] ]
        ch_vcfs = BCFTOOLS_SORT_ME.out.vcf
            .map { meta, vcf ->
                [ groupKey( meta - meta.subMap('interval'), meta.nr_of_intervals ), vcf ]
            }
            .groupTuple()
            .map { meta, vcf ->
                [ meta - meta.subMap('nr_of_intervals'), vcf ]
            }

        ch_tbis = TABIX_ME_SPLIT.out.index
            .map { meta, tbi ->
                [ groupKey( meta - meta.subMap('interval'), meta.nr_of_intervals ), tbi ]
            }
            .groupTuple()
            .map { meta, tbi ->
                [ meta - meta.subMap('nr_of_intervals'), tbi ]
            }

        ch_vcfs_tbis = ch_vcfs.join( ch_tbis, failOnMismatch: true )

        // Concatenate the chromosome vcfs to sample vcfs
        BCFTOOLS_CONCAT_ME ( ch_vcfs_tbis )

        // Merge sample vcfs to a case vcf
        ch_vcf_list = BCFTOOLS_CONCAT_ME.out.vcf
            .collect{_meta, vcf -> vcf}
            .toList()
            .collect()

        ch_svdb_merge_me_input = ch_case_info
            .combine(ch_vcf_list)

        SVDB_MERGE_ME ( ch_svdb_merge_me_input, [], true )
        TABIX_ME ( SVDB_MERGE_ME.out.vcf )

    emit:
        tbi = TABIX_ME.out.index    // channel: [ val(meta), path(tbi) ]
        vcf = SVDB_MERGE_ME.out.vcf // channel: [ val(meta), path(vcf) ]
}
