//
// Call SV MT
//

include { MT_DELETION         } from '../../../modules/local/mt_deletion_script'
include { PREP_MITOSALT       } from '../../../modules/local/prep_mitosalt/main'
include { MITOSALT            } from '../../../modules/local/mitosalt/main'
include { SEQTK_SAMPLE        } from '../../../modules/nf-core/seqtk/sample/main'
include { SALTSHAKER_CALL     } from '../../../modules/nf-core/saltshaker/call/main'
include { SALTSHAKER_CLASSIFY } from '../../../modules/nf-core/saltshaker/classify/main'
include { SALTSHAKER_PLOT     } from '../../../modules/nf-core/saltshaker/plot/main'
include { SVDB_MERGE          } from '../../../modules/nf-core/svdb/merge/main'

workflow CALL_SV_MT {
    take:
        ch_bam_bai                            // channel: [mandatory] [ val(meta), path(bam) ]
        ch_case_info                          // channel: [mandatory] [ val(case_info) ]
        ch_genome_chrsizes                    // channel: [mandatory] [ path(chrsizes) ]
        ch_genome_fai                         // channel: [mandatory] [ val(meta), path(genomefai) ]
        ch_genome_fasta                       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_hisat2index                 // channel: [mandatory] [ val(meta), path(hisat2index) ]
        ch_mt_fai                             // channel: [mandatory] [ val(meta), path(mtfai) ]
        ch_mt_fasta                           // channel: [mandatory] [ val(meta), path(mtfasta) ]
        ch_mt_lastdb                          // channel: [mandatory] [ val(meta), path(lastindex) ]
        ch_reads                              // channel: [mandatory] [ val(meta), [path(reads)] ]
        ch_subdepth                           // channel: [mandatory] [ val(mitosalt_depth) ]
        ch_svcaller_priority                  // channel: [mandatory] [ val(["var caller tag 1", ...]) ]
        val_heavy_strand_origin_start         // string: [mandatory] mitochondira_heavy_strand_origin_start
        val_heavy_strand_origin_end           // string: [mandatory] mitochondira_heavy_strand_origin_end
        val_light_strand_origin_start         // string: [mandatory] mitochondira_light_strand_origin_start
        val_light_strand_origin_end           // string: [mandatory] mitochondira_light_strand_origin_end
        val_mito_length                       // string: [mandatory] mito_length
        val_mito_name                         // string: [mandatory] mito_name
        val_mitosalt_breakspan                // string: [mandatory] mitosalt_breakspan
        val_mitosalt_breakthreshold           // string: [mandatory] mitosalt_breakthreshold
        val_mitosalt_cluster_threshold        // string: [mandatory] mitosalt_cluster_threshold
        val_mitosalt_deletion_threshold_max   // string: [mandatory] mitosalt_deletion_threshold_max
        val_mitosalt_deletion_threshold_min   // string: [mandatory] mitosalt_deletion_threshold_min
        val_mitosalt_evalue_threshold         // string: [mandatory] mitosalt_evalue_threshold
        val_mitosalt_exclude                  // string: [mandatory] mitosalt_exclude
        val_mitosalt_flank                    // string: [mandatory] mitosalt_flank
        val_mitosalt_heteroplasmy_limit       // string: [mandatory] mitosalt_heteroplasmy_limit
        val_mitosalt_paired_distance          // string: [mandatory] mitosalt_paired_distance
        val_mitosalt_score_threshold          // string: [mandatory] mitosalt_score_threshold
        val_mitosalt_sizelimit                // string: [mandatory] mitosalt_sizelimit
        val_mitosalt_split_distance_threshold // string: [mandatory] mitosalt_split_distance_threshold
        val_mitosalt_split_length             // string: [mandatory] mitosalt_split_length

    main:
        ch_saltshaker_txt   = channel.empty()
        ch_saltshaker_vcf   = channel.empty()
        ch_saltshaker_plot  = channel.empty()

        if (!(params.skip_tools && params.skip_tools.split(',').contains('mitosalt'))) {
            ch_reads_subdepth      = ch_reads.combine(ch_subdepth)

            SEQTK_SAMPLE (ch_reads_subdepth)

            PREP_MITOSALT(
                ch_genome_chrsizes,
                ch_genome_fai,
                ch_genome_hisat2index,
                ch_mt_fai,
                ch_mt_fasta,
                ch_mt_lastdb,
                val_mitosalt_breakspan,
                val_mitosalt_breakthreshold,
                val_mitosalt_cluster_threshold,
                val_mitosalt_deletion_threshold_max,
                val_mitosalt_deletion_threshold_min,
                val_mitosalt_evalue_threshold,
                val_mitosalt_exclude,
                val_mitosalt_flank,
                val_mitosalt_heteroplasmy_limit,
                val_mito_name,
                val_mitosalt_paired_distance,
                val_mitosalt_score_threshold,
                val_mitosalt_sizelimit,
                val_mitosalt_split_distance_threshold,
                val_mitosalt_split_length
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

            MITOSALT.out.cluster
                .filter{ _meta, out -> out.countLines() > 0 }
                .set{ch_cluster}

            MITOSALT.out.breakpoint
                .join(ch_cluster)
                .set{ch_saltshaker_in}

            SALTSHAKER_CALL(
                ch_saltshaker_in,
                ch_mt_fasta,
                val_mitosalt_flank,
                val_mitosalt_heteroplasmy_limit,
                val_mito_length,
                val_heavy_strand_origin_start,
                val_heavy_strand_origin_end,
                val_light_strand_origin_start,
                val_light_strand_origin_end
            )

            SALTSHAKER_CLASSIFY(
                SALTSHAKER_CALL.out.call,
                val_mito_name
            )
            ch_saltshaker_txt = SALTSHAKER_CLASSIFY.out.txt
            ch_saltshaker_vcf = SALTSHAKER_CLASSIFY.out.vcf

            SALTSHAKER_PLOT(
                SALTSHAKER_CLASSIFY.out.classify
            )
            ch_saltshaker_plot = SALTSHAKER_PLOT.out.plot

            SALTSHAKER_CLASSIFY.out.vcf
                .collect{_meta, vcf -> vcf}
                .toList()
                .set { ch_vcf_file_list }

            ch_case_info
                .combine(ch_vcf_file_list)
                .filter{ it -> it.size() == 2}
                .set { ch_merge_input_vcfs }

            SVDB_MERGE ( ch_merge_input_vcfs, [], true ).vcf
                .set {ch_saltshaker_vcf}
                // Update priority list when we know saltshaker will run
            ch_svcaller_priority = ch_svcaller_priority
                .concat(ch_saltshaker_vcf.map{ _ -> ["mitosalt"] })
                .collect()

        }
        MT_DELETION(ch_bam_bai, ch_genome_fasta)

        ch_publish = ch_saltshaker_txt
                        .mix(ch_saltshaker_vcf)
                        .mix(ch_saltshaker_plot)
                        .mix(MT_DELETION.out.mt_del_result)
                        .map { meta, value -> ['call_sv/', [meta, value]] }

    emit:
        saltshaker_txt   = ch_saltshaker_txt             // channel: [ val(meta), path(txt) ]
        saltshaker_vcf   = ch_saltshaker_vcf             // channel: [ val(meta), path(vcf) ]
        saltshaker_plot  = ch_saltshaker_plot            // channel: [ val(meta), path(png) ]
        mt_del_result    = MT_DELETION.out.mt_del_result // channel: [ val(meta), path(txt) ]
        updated_priority = ch_svcaller_priority          // channel: [ val(["caller1", "caller2", ...]) ] - includes "mitosalt" if it ran
        publish          = ch_publish                    // channel: [ val(destination), val(value) ]
}
