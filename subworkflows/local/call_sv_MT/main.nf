//
// Call SV MT
//

include { MT_DELETION   } from '../../../modules/local/mt_deletion_script'
include { PREP_MITOSALT } from '../../../modules/local/prep_mitosalt/main'
include { MITOSALT      } from '../../../modules/local/mitosalt/main'
include { SEQTK_SAMPLE  } from '../../../modules/nf-core/seqtk/sample/main'
include { CAT_FASTQ     } from '../../../modules/nf-core/cat/fastq/main'

workflow CALL_SV_MT {
    take:
        ch_bam_bai                   // channel: [mandatory] [ val(meta), path(bam) ]
        ch_genome_chrsizes           // channel: [mandatory] [ path(chrsizes) ]
        ch_genome_fai                // channel: [mandatory] [ val(meta), path(genomefai) ]
        ch_genome_fasta              // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_hisat2index        // channel: [mandatory] [ val(meta), path(hisat2index) ]
        ch_mt_fai                    // channel: [mandatory] [ val(meta), path(mtfai) ]
        ch_mt_fasta                  // channel: [mandatory] [ val(meta), path(mtfasta) ]
        ch_mt_lastdb                 // channel: [mandatory] [ val(meta), path(lastindex) ]
        ch_reads                     // channel: [mandatory] [ val(meta), [path(reads)] ]
        ch_subdepth                  // channel: [mandatory] [ val(mitosalt_depth) ]
        val_breakspan                // string: [mandatory] mitosalt_break_span
        val_breakthreshold           // string: [mandatory] mitosalt_break_threshold
        val_cluster_threshold        // string: [mandatory] mitosalt_cluster_threshold
        val_deletion_threshold_max   // string: [mandatory] mitosalt_del_threshold_max
        val_deletion_threshold_min   // string: [mandatory] mitosalt_del_threshold_min
        val_evalue_threshold         // string: [mandatory] mitosalt_evalue_threshold
        val_exclude                  // string: [mandatory] mitosalt_exclude
        val_flank                    // string: [mandatory] mitosalt_flank
        val_hplimit                  // string: [mandatory] mitosalt_hp_limit
        val_mito_name                // string: [mandatory] mito_name
        val_paired_distance          // string: [mandatory] mitosalt_paired_distance
        val_score_threshold          // string: [mandatory] mitosalt_score_threshold
        val_sizelimit                // string: [mandatory] mitosalt_size_limit
        val_split_distance_threshold // string: [mandatory] mitosalt_split_dist_threshold
        val_split_length             // string: [mandatory] mitosalt_split_length

    main:
        ch_mitosalt_publish = channel.empty()

        if (!(params.skip_tools && params.skip_tools.split(',').contains('mitosalt'))) {

            ch_reads
                .map { meta, reads ->
                        def groupKey = meta.sample
                        return [groupKey, meta, reads]
                }
                .groupTuple(by: 0)
                .map { sample_id, meta_list, reads_list ->
                    def combined_meta = meta_list[0].clone()
                    combined_meta.id = sample_id
                    combined_meta.remove('lane')
                    combined_meta.remove('read_group')

                    def all_reads = reads_list.flatten()

                    [combined_meta, all_reads]
                }
                .set {ch_cat_fastq}

            CAT_FASTQ(ch_cat_fastq)

            ch_reads_subdepth = CAT_FASTQ.out.reads.combine(ch_subdepth)

            SEQTK_SAMPLE (ch_reads_subdepth)

            PREP_MITOSALT(
                ch_genome_chrsizes,
                ch_genome_fai,
                ch_genome_hisat2index,
                ch_mt_fai,
                ch_mt_fasta,
                ch_mt_lastdb,
                val_breakspan,
                val_breakthreshold,
                val_cluster_threshold,
                val_deletion_threshold_max,
                val_deletion_threshold_min,
                val_evalue_threshold,
                val_exclude,
                val_flank,
                val_hplimit,
                val_mito_name,
                val_paired_distance,
                val_score_threshold,
                val_sizelimit,
                val_split_distance_threshold,
                val_split_length
            )

            MITOSALT(
                SEQTK_SAMPLE.out.reads,
                PREP_MITOSALT.out.msconfig,
                ch_genome_chrsizes,
                ch_genome_fai,
                ch_genome_hisat2index,
                ch_mt_fai,
                ch_mt_fasta,
                ch_mt_lastdb
            )
            ch_mitosalt_publish = MITOSALT.out.breakpoint
                .mix(MITOSALT.out.cluster)

        }
        MT_DELETION(ch_bam_bai, ch_genome_fasta)

        ch_publish = ch_mitosalt_publish
            .mix(MT_DELETION.out.mt_del_result)
            .map { meta, value -> ['call_sv/mitochondria/', [meta, value]] }

    emit:
        mitosalt_breakpoint = MITOSALT.out.breakpoint       // channel: [ val(meta), path(breakpoint) ]
        mitosalt_cluster    = MITOSALT.out.cluster          // channel: [ val(meta), path(cluster) ]
        mt_del_result       = MT_DELETION.out.mt_del_result // channel: [ val(meta), path(txt) ]
        publish = ch_publish                                // channel: [ val(destination), val(value) ]
}
