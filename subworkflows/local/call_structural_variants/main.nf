//
// A nested subworkflow to call structural variants.
//

include { CALL_SV_MANTA                  } from '../call_sv_manta'
include { CALL_SV_TIDDIT                 } from '../call_sv_tiddit'
include { SVDB_MERGE                     } from '../../../modules/nf-core/svdb/merge/main'
include { CALL_SV_GERMLINECNVCALLER      } from '../call_sv_germlinecnvcaller'
include { CALL_SV_CNVNATOR               } from '../call_sv_cnvnator'
include { TABIX_TABIX                    } from '../../../modules/nf-core/tabix/tabix/main'

workflow CALL_STRUCTURAL_VARIANTS {

    take:
        ch_genome_bam          // channel: [mandatory] [ val(meta), path(bam) ]
        ch_genome_bai          // channel: [mandatory] [ val(meta), path(bai) ]
        ch_genome_bam_bai      // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_bwa_index           // channel: [mandatory] [ val(meta), path(index)]
        ch_genome_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai          // channel: [mandatory] [ val(meta), path(fai) ]
        ch_case_info           // channel: [mandatory] [ val(case_info) ]
        ch_target_bed          // channel: [mandatory for WES] [ val(meta), path(bed), path(tbi) ]
        ch_genome_dictionary   // channel: [optional; used by mandatory for GATK's cnvcaller][ val(meta), path(dict) ]
        ch_svcaller_priority   // channel: [mandatory] [ val(["var caller tag 1", ...]) ]
        ch_readcount_intervals // channel: [optional; used by mandatory for GATK's cnvcaller][ path(intervals) ]
        ch_ploidy_model        // channel: [optional; used by mandatory for GATK's cnvcaller][ path(ploidy_model) ]
        ch_gcnvcaller_model    // channel: [optional; used by mandatory for GATK's cnvcaller][ path(gcnvcaller_model) ]
        val_analysis_type      // string: "wes", "wgs", or "mito"
        skip_germlinecnvcaller // boolean

    main:
        ch_versions   = channel.empty()
        ch_merged_svs = channel.empty()
        ch_merged_tbi = channel.empty()

        if (!val_analysis_type.equals("mito")) {
            CALL_SV_MANTA (ch_genome_bam, ch_genome_bai, ch_genome_fasta, ch_genome_fai, ch_case_info, ch_target_bed, val_analysis_type)
                .filtered_diploid_sv_vcf
                .collect{ _meta, vcf -> vcf }
                .set{ manta_vcf }
        }

        if (val_analysis_type.equals("wgs")) {
            CALL_SV_TIDDIT (ch_genome_bam_bai, ch_genome_fasta, ch_bwa_index, ch_case_info)
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

            ch_versions = ch_versions.mix(CALL_SV_GERMLINECNVCALLER.out.versions)
        }

        //merge
        if (skip_germlinecnvcaller) {
            if (val_analysis_type.equals("wgs")) {
                tiddit_vcf
                    .combine(manta_vcf)
                    .combine(cnvnator_vcf)
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
                .toList()
                .set { vcf_list }
        } else if (!val_analysis_type.equals("mito")) {
            manta_vcf
                .combine(gcnvcaller_vcf)
                .toList()
                .set { vcf_list }
        }

        if (!val_analysis_type.equals("mito")) {
            ch_case_info
                .combine(vcf_list)
                .set { merge_input_vcfs }

            SVDB_MERGE (merge_input_vcfs, ch_svcaller_priority, true)

            TABIX_TABIX (SVDB_MERGE.out.vcf)
            ch_merged_svs = SVDB_MERGE.out.vcf
            ch_merged_tbi = TABIX_TABIX.out.index
        }

    emit:
        vcf      = ch_merged_svs // channel: [ val(meta), path(vcf)]
        tbi      = ch_merged_tbi // channel: [ val(meta), path(tbi)]
        versions = ch_versions   // channel: [ path(versions.yml) ]
}
