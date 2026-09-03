//
// A nested subworkflow to call structural variants from nuclear DNA.
//

include { CALL_SV_CNVNATOR          } from '../call_sv_cnvnator'
include { CALL_SV_GERMLINECNVCALLER } from '../call_sv_germlinecnvcaller'
include { CALL_SV_MANTA             } from '../call_sv_manta'
include { CALL_SV_TIDDIT            } from '../call_sv_tiddit'

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
        skip_germlinecnvcaller                // boolean
        val_analysis_type                     // string: "wes", "wgs", or "mito"

    main:
        ch_cnvnator_vcf    = channel.empty()
        ch_gcnvcaller_vcf  = channel.empty()
        ch_manta_vcf       = channel.empty()
        ch_tiddit_vcf      = channel.empty()

        // CALL_SV is only invoked for non-mito analysis types (gated at the call site in
        // raredisease.nf, mirroring CALL_SV_MT's val_run_mt gate), so no mito check is needed here.
        ch_manta_vcf = CALL_SV_MANTA (ch_genome_bam, ch_genome_bai, ch_genome_fasta, ch_genome_fai, ch_case_info, ch_manta_regions)
            .filtered_diploid_sv_vcf
            .collect{ _meta, vcf -> vcf }

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

        // Collect individual caller VCFs in a fixed, consistent order so the flat list lines
        // up positionally with the svcaller_priority tags built in main.nf: tiddit -> manta ->
        // gcnvcaller -> cnvnator.
        // Empty channels won't contribute any items.
        ch_vcf_paths = ch_tiddit_vcf
            .concat(ch_manta_vcf)
            .concat(ch_gcnvcaller_vcf)
            .concat(ch_cnvnator_vcf)
            .collect()

    emit:
        vcfs = ch_vcf_paths // channel: [ [path(vcf), path(vcf), ...] ] - flat list of nuclear caller VCFs, in priority order
}
