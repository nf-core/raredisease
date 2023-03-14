//
// A nested subworkflow to call structural variants.
//

include { CALL_SV_MANTA     } from './variant_calling/call_sv_manta'
include { CALL_SV_TIDDIT    } from './variant_calling/call_sv_tiddit'
include { SVDB_MERGE        } from '../../modules/nf-core/svdb/merge/main'
include { CALL_CNV_CNVPYTOR } from './variant_calling/call_cnv_cnvpytor'

workflow CALL_STRUCTURAL_VARIANTS {

    take:
        ch_bam            // channel: [mandatory] [ val(meta), path(bam) ]
        ch_bai            // channel: [mandatory] [ val(meta), path(bai) ]
        ch_bam_bai        // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_bwa_index      // channel: [mandatory] [ val(meta), path(index)]
        ch_fasta_no_meta  // channel: [mandatory] [ path(fasta) ]
        ch_fasta_meta     // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_fai            // channel: [mandatory] [ path(fai) ]
        ch_case_info      // channel: [mandatory] [ val(case_info) ]
        ch_target_bed     // channel: [mandatory for WES] [ val(meta), path(bed), path(tbi) ]
        val_cnvpytor_bins // string: [optional] binsizes for cnvpytor default: 1000

    main:
        ch_versions = Channel.empty()

        CALL_SV_MANTA (ch_bam, ch_bai, ch_fasta_no_meta, ch_fai, ch_case_info, ch_target_bed)
            .diploid_sv_vcf
            .collect{it[1]}
            .set{ manta_vcf }

        CALL_SV_TIDDIT (ch_bam_bai, ch_fasta_meta, ch_bwa_index, ch_case_info)
            .vcf
            .collect{it[1]}
            .set { tiddit_vcf }

        //cnvpytor
        // CALL_CNV_CNVPYTOR ( bam, bai, case_info, cnvpytor_bins, fasta_no_meta, fai)
        //     .candidate_cnvs_vcf
        //     .collect{it[1]}
        //     .set {cnvpytor_vcf }

        // //merge
        // tiddit_vcf
        //     .combine(manta_vcf)
        //     .combine(cnvpytor_vcf)
        //     .toList()
        //     .set { vcf_list }

        //merge
        tiddit_vcf
            .combine(manta_vcf)
            .toList()
            .set { vcf_list }

        ch_case_info
            .combine(vcf_list)
            .set { merge_input_vcfs }

        // SVDB_MERGE ( merge_input_vcfs, ["tiddit","manta","cnvpytor"] )
        SVDB_MERGE (merge_input_vcfs, ["tiddit","manta"])

        ch_versions = ch_versions.mix(CALL_SV_MANTA.out.versions)
        ch_versions = ch_versions.mix(CALL_SV_TIDDIT.out.versions)
        // ch_versions = ch_versions.mix(CALL_CNV_CNVPYTOR.out.versions)

    emit:
        vcf      = SVDB_MERGE.out.vcf // channel: [ val(meta), path(vcf)]
        versions = ch_versions        // channel: [ path(versions.yml) ]
}
