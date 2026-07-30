//
// A nested subworkflow to call structural variants from nuclear DNA.
//

include { CALL_SV_CNVNATOR          } from '../call_sv_cnvnator'
include { CALL_SV_GERMLINECNVCALLER } from '../call_sv_germlinecnvcaller'
include { CALL_SV_MANTA             } from '../call_sv_manta'
include { CALL_SV_TIDDIT            } from '../call_sv_tiddit'
include { SVDB_MERGE                } from '../../../modules/nf-core/svdb/merge/main'
include { TABIX_TABIX               } from '../../../modules/nf-core/tabix/tabix/main'

workflow CALL_SV {

    take:
        ch_bwa_index                          // channel: [mandatory] [ val(meta), path(index)]
        ch_case_info                          // channel: [mandatory] [ val(case_info) ]
        ch_gcnvcaller_model                   // channel: [optional; used by mandatory for GATK's cnvcaller][ path(gcnvcaller_model) ]
        ch_genome_bai                         // channel: [mandatory] [ val(meta), path(bai) ]
        ch_genome_bam                         // channel: [mandatory] [ val(meta), path(bam) ]
        ch_genome_bam_bai                     // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_dictionary                  // channel: [optional; used by mandatory for GATK's cnvcaller][ val(meta), path(dict) ]
        ch_genome_fai                         // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_fasta                       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_manta_regions                      // channel: [mandatory] [ path(bed), path(tbi) ]
        ch_ploidy_model                       // channel: [optional; used by mandatory for GATK's cnvcaller][ path(ploidy_model) ]
        ch_readcount_intervals                // channel: [optional; used by mandatory for GATK's cnvcaller][ path(intervals) ]
        ch_svcaller_priority                  // channel: [mandatory] [ val(["var caller tag 1", ...]) ]
        skip_germlinecnvcaller                // boolean
        val_analysis_type                     // string: "wes", "wgs", or "mito"

    main:
        ch_cnvnator_vcf    = channel.empty()
        ch_gcnvcaller_vcf  = channel.empty()
        ch_manta_vcf       = channel.empty()
        ch_merged_svs      = channel.empty()
        ch_merged_tbi      = channel.empty()
        ch_tiddit_vcf      = channel.empty()

        if (!val_analysis_type.equals("mito")) {
            ch_manta_vcf = CALL_SV_MANTA (ch_genome_bam, ch_genome_bai, ch_genome_fasta, ch_genome_fai, ch_case_info, ch_manta_regions)
                .filtered_diploid_sv_vcf
                .collect{ _meta, vcf -> vcf }
        }

        if (val_analysis_type.equals("wgs")) {
            ch_tiddit_vcf = CALL_SV_TIDDIT (ch_genome_bam_bai, ch_genome_fai, ch_genome_fasta, ch_bwa_index, ch_case_info)
                .vcf
                .collect{ _meta, vcf -> vcf }

            ch_cnvnator_vcf = CALL_SV_CNVNATOR (ch_genome_bam_bai, ch_genome_fasta, ch_case_info)
                .vcf
                .collect{ _meta, vcf -> vcf }
        }

        if (!skip_germlinecnvcaller) {
            ch_gcnvcaller_vcf = CALL_SV_GERMLINECNVCALLER (ch_genome_bam_bai, ch_genome_fasta, ch_genome_fai, ch_readcount_intervals, ch_genome_dictionary, ch_ploidy_model, ch_gcnvcaller_model, ch_case_info)
                .genotyped_filtered_segments_vcf
                .collect{ _meta, vcf -> vcf }

        }

        // Merge - with consistent ordering using concat. Only meaningful outside of mito-only
        // analysis, since none of the nuclear callers above run for that mode.
        if (!val_analysis_type.equals("mito")) {
            // Concatenate in specific order: tiddit -> manta -> gcnvcaller -> cnvnator
            // Empty channels won't contribute any items
            ch_vcf_paths = ch_tiddit_vcf
                .concat(ch_manta_vcf)
                .concat(ch_gcnvcaller_vcf)
                .concat(ch_cnvnator_vcf)
                .collect()
                .map { vcf_list -> [vcf_list] }
            ch_merge_vcfs_in = ch_case_info
                .combine(ch_vcf_paths)
            SVDB_MERGE (ch_merge_vcfs_in, ch_svcaller_priority, false)

            TABIX_TABIX (SVDB_MERGE.out.vcf)
            ch_merged_svs = SVDB_MERGE.out.vcf
            ch_merged_tbi = TABIX_TABIX.out.index
        }

    emit:
        tbi = ch_merged_tbi // channel: [ val(meta), path(tbi)]
        vcf = ch_merged_svs // channel: [ val(meta), path(vcf)]
}
