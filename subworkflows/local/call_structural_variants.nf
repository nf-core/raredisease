//
// A nested subworkflow to call structural variants.
//

include { CALL_SV_MANTA                  } from './variant_calling/call_sv_manta'
include { CALL_SV_MT                     } from './variant_calling/call_sv_MT'
include { CALL_SV_MT as CALL_SV_MT_SHIFT } from './variant_calling/call_sv_MT'
include { CALL_SV_TIDDIT                 } from './variant_calling/call_sv_tiddit'
include { SVDB_MERGE                     } from '../../modules/nf-core/svdb/merge/main'
include { CALL_SV_GERMLINECNVCALLER      } from './variant_calling/call_sv_germlinecnvcaller'
include { CALL_SV_CNVNATOR               } from './variant_calling/call_sv_cnvnator'
include { TABIX_TABIX                    } from '../../modules/nf-core/tabix/tabix/main'

workflow CALL_STRUCTURAL_VARIANTS {

    take:
        ch_genome_bam          // channel: [mandatory] [ val(meta), path(bam) ]
        ch_genome_bai          // channel: [mandatory] [ val(meta), path(bai) ]
        ch_genome_bam_bai      // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_mt_bam_bai          // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_mtshift_bam_bai     // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_bwa_index           // channel: [mandatory] [ val(meta), path(index)]
        ch_genome_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai          // channel: [mandatory] [ val(meta), path(fai) ]
        ch_mtshift_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_case_info           // channel: [mandatory] [ val(case_info) ]
        ch_target_bed          // channel: [mandatory for WES] [ val(meta), path(bed), path(tbi) ]
        ch_genome_dictionary   // channel: [optional; used by mandatory for GATK's cnvcaller][ val(meta), path(dict) ]
        ch_svcaller_priority   // channel: [mandatory] [ val(["var caller tag 1", ...]) ]
        ch_readcount_intervals // channel: [optional; used by mandatory for GATK's cnvcaller][ path(intervals) ]
        ch_ploidy_model        // channel: [optional; used by mandatory for GATK's cnvcaller][ path(ploidy_model) ]
        ch_gcnvcaller_model    // channel: [optional; used by mandatory for GATK's cnvcaller][ path(gcnvcaller_model) ]

    main:
        ch_versions = Channel.empty()

        CALL_SV_MANTA (ch_genome_bam, ch_genome_bai, ch_genome_fasta, ch_genome_fai, ch_case_info, ch_target_bed)
            .diploid_sv_vcf
            .collect{it[1]}
            .set{ manta_vcf }

        CALL_SV_TIDDIT (ch_genome_bam_bai, ch_genome_fasta, ch_bwa_index, ch_case_info)
            .vcf
            .collect{it[1]}
            .set { tiddit_vcf }

        if (!params.skip_germlinecnvcaller) {
            CALL_SV_GERMLINECNVCALLER (ch_genome_bam_bai, ch_genome_fasta, ch_genome_fai, ch_readcount_intervals, ch_genome_dictionary, ch_ploidy_model, ch_gcnvcaller_model)
                .genotyped_filtered_segments_vcf
                .collect{it[1]}
                .set { gcnvcaller_vcf }

            ch_versions = ch_versions.mix(CALL_SV_GERMLINECNVCALLER.out.versions)
        }

        CALL_SV_CNVNATOR (ch_genome_bam_bai, ch_genome_fasta, ch_genome_fai, ch_case_info)
            .vcf
            .collect{it[1]}
            .set { cnvnator_vcf }

        CALL_SV_MT (ch_mt_bam_bai, ch_genome_fasta)

        //merge
        if (params.skip_germlinecnvcaller) {
            tiddit_vcf
                .combine(manta_vcf)
                .combine(cnvnator_vcf)
                .toList()
                .set { vcf_list }
        } else {
            tiddit_vcf
                .combine(manta_vcf)
                .combine(gcnvcaller_vcf)
                .combine(cnvnator_vcf)
                .toList()
                .set { vcf_list }
        }

        ch_case_info
            .combine(vcf_list)
            .set { merge_input_vcfs }

        SVDB_MERGE (merge_input_vcfs, ch_svcaller_priority)

        TABIX_TABIX (SVDB_MERGE.out.vcf)

        ch_versions = ch_versions.mix(CALL_SV_MANTA.out.versions)
        ch_versions = ch_versions.mix(CALL_SV_MT.out.versions)
        ch_versions = ch_versions.mix(CALL_SV_TIDDIT.out.versions)
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
        vcf      = SVDB_MERGE.out.vcf  // channel: [ val(meta), path(vcf)]
        tbi      = TABIX_TABIX.out.tbi // channel: [ val(meta), path(tbi)]
        versions = ch_versions         // channel: [ path(versions.yml) ]
}
