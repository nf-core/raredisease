//
// A nested subworkflow to call structural variants.
//

include { CALL_SV_MANTA                  } from '../call_sv_manta'
include { CALL_SV_TIDDIT                 } from '../call_sv_tiddit'
include { CALL_SV_MT                     } from '../call_sv_MT'
include { SVDB_MERGE                     } from '../../../modules/nf-core/svdb/merge/main'
include { CALL_SV_GERMLINECNVCALLER      } from '../call_sv_germlinecnvcaller'
include { CALL_SV_CNVNATOR               } from '../call_sv_cnvnator'
include { TABIX_TABIX                    } from '../../../modules/nf-core/tabix/tabix/main'

workflow CALL_STRUCTURAL_VARIANTS {

    take:
        ch_genome_bam                         // channel: [mandatory] [ val(meta), path(bam) ]
        ch_genome_bai                         // channel: [mandatory] [ val(meta), path(bai) ]
        ch_genome_bam_bai                     // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_bwa_index                          // channel: [mandatory] [ val(meta), path(index)]
        ch_genome_fasta                       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai                         // channel: [mandatory] [ val(meta), path(fai) ]
        ch_case_info                          // channel: [mandatory] [ val(case_info) ]
        ch_target_bed                         // channel: [mandatory for WES] [ val(meta), path(bed), path(tbi) ]
        ch_genome_dictionary                  // channel: [optional; used by mandatory for GATK's cnvcaller][ val(meta), path(dict) ]
        ch_svcaller_priority                  // channel: [mandatory] [ val(["var caller tag 1", ...]) ]
        ch_readcount_intervals                // channel: [optional; used by mandatory for GATK's cnvcaller][ path(intervals) ]
        ch_ploidy_model                       // channel: [optional; used by mandatory for GATK's cnvcaller][ path(ploidy_model) ]
        ch_gcnvcaller_model                   // channel: [optional; used by mandatory for GATK's cnvcaller][ path(gcnvcaller_model) ]
        val_analysis_type                     // string: "wes", "wgs", or "mito"
        skip_germlinecnvcaller                // boolean
        ch_mt_bam_bai                         // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_chrsizes                    // channel: [mandatory] [ path(chrsizes) ]
        ch_genome_hisat2index                 // channel: [mandatory] [ val(meta), path(hisat2index) ]
        ch_mt_fai                             // channel: [mandatory] [ val(meta), path(mtfai) ]
        ch_mt_fasta                           // channel: [mandatory] [ val(meta), path(mtfasta) ]
        ch_mt_lastdb                          // channel: [mandatory] [ val(meta), path(lastindex) ]
        ch_reads                              // channel: [mandatory] [ val(meta), [path(reads)] ]
        ch_subdepth                           // channel: [mandatory] [ val(mitosalt_depth) ]
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
        val_run_mt_for_wes                    // boolean: [mandatory] run_mt_for_wes

    main:
        ch_merged_svs = channel.empty()
        ch_merged_tbi = channel.empty()

        if (!val_analysis_type.equals("mito")) {
            CALL_SV_MANTA (ch_genome_bam, ch_genome_bai, ch_genome_fasta, ch_genome_fai, ch_case_info, ch_target_bed, val_analysis_type)
                .filtered_diploid_sv_vcf
                .collect{ _meta, vcf -> vcf }
                .set{ manta_vcf }
        }

        if (val_analysis_type.equals("wgs")) {
            CALL_SV_TIDDIT (ch_genome_bam_bai, ch_genome_fai, ch_genome_fasta, ch_bwa_index, ch_case_info)
                .vcf
                .collect{ _meta, vcf -> vcf }
                .set { tiddit_vcf }

            CALL_SV_CNVNATOR (ch_genome_bam_bai, ch_genome_fasta, ch_genome_fai, ch_case_info)
                .vcf
                .collect{ _meta, vcf -> vcf }
                .set { cnvnator_vcf }
        }

        if (!skip_germlinecnvcaller) {
            CALL_SV_GERMLINECNVCALLER (ch_genome_bam_bai, ch_genome_fasta, ch_genome_fai, ch_readcount_intervals, ch_genome_dictionary, ch_ploidy_model, ch_gcnvcaller_model, ch_case_info)
                .genotyped_filtered_segments_vcf
                .collect{ _meta, vcf -> vcf }
                .set { gcnvcaller_vcf }

        }

        if (val_analysis_type.matches("wgs|mito") || val_run_mt_for_wes) {
            CALL_SV_MT(
                ch_mt_bam_bai,
                ch_genome_chrsizes,
                ch_genome_fai,
                ch_genome_fasta,
                ch_genome_hisat2index,
                ch_mt_fai,
                ch_mt_fasta,
                ch_mt_lastdb,
                ch_reads,
                ch_subdepth,
                val_heavy_strand_origin_start,
                val_heavy_strand_origin_end,
                val_light_strand_origin_start,
                val_light_strand_origin_end,
                val_mito_length,
                val_mito_name,
                val_mitosalt_breakspan,
                val_mitosalt_breakthreshold,
                val_mitosalt_cluster_threshold,
                val_mitosalt_deletion_threshold_max,
                val_mitosalt_deletion_threshold_min,
                val_mitosalt_evalue_threshold,
                val_mitosalt_exclude,
                val_mitosalt_flank,
                val_mitosalt_heteroplasmy_limit,
                val_mitosalt_paired_distance,
                val_mitosalt_score_threshold,
                val_mitosalt_sizelimit,
                val_mitosalt_split_distance_threshold,
                val_mitosalt_split_length)
                .mitosalt_vcf
                .collect{ _meta, vcf -> vcf }
                .set { mitosalt_vcf }
        }

        //merge
        if (skip_germlinecnvcaller) {
            if (val_analysis_type.equals("wgs")) {
                tiddit_vcf
                    .combine(manta_vcf)
                    .combine(cnvnator_vcf)
                    .combine(mitosalt_vcf)
                    .toList()
                    .set { vcf_list }
            } else if (!val_analysis_type.equals("mito")) {
                manta_vcf
                    .toList()
                    .set { vcf_list }
            }
        } else if (val_analysis_type.equals("wgs")) {
            tiddit_vcf
                .combine(manta_vcf)
                .combine(gcnvcaller_vcf)
                .combine(cnvnator_vcf)
                .combine(mitosalt_vcf)
                .toList()
                .set { vcf_list }
        } else if (!val_analysis_type.equals("mito")) {
            manta_vcf
                .combine(gcnvcaller_vcf)
                .toList()
                .set { vcf_list }
        }

        if (!val_analysis_type.equals("mito")) {
            ch_case_info.view()
            ch_svcaller_priority.view()
            vcf_list.view()
            ch_case_info
                .combine(vcf_list)
                .set { merge_input_vcfs }

            SVDB_MERGE (merge_input_vcfs, ch_svcaller_priority, true)

            TABIX_TABIX (SVDB_MERGE.out.vcf)
            ch_merged_svs = SVDB_MERGE.out.vcf
            ch_merged_tbi = TABIX_TABIX.out.index
        } else {
            TABIX_TABIX (mitosalt_vcf)
            ch_merged_svs = mitosalt_vcf
            ch_merged_tbi = TABIX_TABIX.out.index
        }

        ch_publish = ch_merged_svs
            .mix(ch_merged_tbi)
            .map { meta, value -> ['call_sv/genome/', [meta, value]] }

    emit:
        vcf      = ch_merged_svs // channel: [ val(meta), path(vcf)]
        tbi      = ch_merged_tbi // channel: [ val(meta), path(tbi)]
        publish = ch_publish     // channel: [ val(destination), val(value) ]
}
