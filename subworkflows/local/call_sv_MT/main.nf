//
// Call SV MT
//

include { MT_DELETION         } from '../../../modules/local/mt_deletion_script'
include { PREP_MITOSALT       } from '../../../modules/local/prep_mitosalt/main'
include { MITOSALT            } from '../../../modules/local/mitosalt/main'
include { SEQTK_SAMPLE        } from '../../../modules/nf-core/seqtk/sample/main'
include { SALTSHAKER_CALL     } from '../../../modules/local/saltshaker/call/main'
include { SALTSHAKER_CLASSIFY } from '../../../modules/local/saltshaker/classify/main'
include { SALTSHAKER_PLOT     } from '../../../modules/local/saltshaker/plot/main'

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
        val_dom_frac                 // string: [mandatory] saltshaker_dominant_fraction
        val_evalue_threshold         // string: [mandatory] mitosalt_evalue_threshold
        val_exclude                  // string: [mandatory] mitosalt_exclude
        val_flank                    // string: [mandatory] mitosalt_flank
        val_group_radius             // string: [mandatory] saltshaker_group_radius
        val_high_het                 // string: [mandatory] saltshaker_high_heteroplasmy
        val_hplimit                  // string: [mandatory] mitosalt_hp_limit
        val_mito_length              // string: [mandatory] mito_length
        val_mito_name                // string: [mandatory] mito_name
        val_mult_thresh              // string: [mandatory] saltshaker_multi_threshold
        val_noise_thresh             // string: [mandatory] saltshaker_noise_threshold
        val_ori_h_start              // string: [mandatory] mito_ori_h_start
        val_ori_h_end                // string: [mandatory] mito_ori_h_end
        val_ori_l_start              // string: [mandatory] mito_ori_l_start
        val_ori_l_end                // string: [mandatory] mito_ori_l_end
        val_paired_distance          // string: [mandatory] mitosalt_paired_distance
        val_score_threshold          // string: [mandatory] mitosalt_score_threshold
        val_sizelimit                // string: [mandatory] mitosalt_size_limit
        val_split_distance_threshold // string: [mandatory] mitosalt_split_dist_threshold
        val_split_length             // string: [mandatory] mitosalt_split_length

    main:
        ch_versions = Channel.empty()

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

            SALTSHAKER_CALL(
                MITOSALT.out.breakpoint,
                MITOSALT.out.cluster,
                ch_mt_fasta,
                val_flank,
                val_hplimit,
                val_mito_length,
                val_ori_h_start,
                val_ori_h_end,
                val_ori_l_start,
                val_ori_l_end
            )

            SALTSHAKER_CLASSIFY(
                SALTSHAKER_CALL.out.call,
                val_dom_frac,
                val_group_radius,
                val_high_het,
                val_mult_thresh,
                val_noise_thresh
            )
            ch_mitosalt_vcf = SALTSHAKER_CLASSIFY.out.classify

            SALTSHAKER_PLOT(
                SALTSHAKER_CLASSIFY.out.classify
            )
            ch_mitosalt_plot = SALTSHAKER_PLOT.out.plot
            

            ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions)
        }
        MT_DELETION(ch_bam_bai, ch_genome_fasta)

    emit:
        mitosalt_vcf        = ch_mitosalt_vcf               // channel: [ val(meta), path(vcf) ]
        mitosalt_plot       = ch_mitosalt_plot              // channel: [ val(meta), path(png) ]
        mt_del_result       = MT_DELETION.out.mt_del_result // channel: [ val(meta), path(txt) ]
        versions            = ch_versions                   // channel: [ path(versions.yml) ]
}
