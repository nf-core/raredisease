//
// A nested subworkflow to call structural variants.
//

include { CALL_SV_MANTA              } from './variant_calling/call_sv_manta'
include { CALL_SV_TIDDIT             } from './variant_calling/call_sv_tiddit'
include { SVDB_MERGE                 } from '../../modules/nf-core/svdb/merge/main'
include { CALL_CNV_GERMLINECNVCALLER } from './variant_calling/call_cnv_germlinecnvcaller'
include { TABIX_TABIX                } from '../../modules/nf-core/tabix/tabix/main'

workflow CALL_STRUCTURAL_VARIANTS {

    take:
        ch_bam            // channel: [mandatory] [ val(meta), path(bam) ]
        ch_bai            // channel: [mandatory] [ val(meta), path(bai) ]
        ch_bam_bai        // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        bam_bai_2     // channel: [ val(meta), path(bam), path(bai) ]
        ch_bwa_index      // channel: [mandatory] [ val(meta), path(index)]
        ch_fasta_no_meta  // channel: [mandatory] [ path(fasta) ]
        ch_fasta_meta     // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_fai            // channel: [mandatory] [ path(fai) ]
        ch_case_info      // channel: [mandatory] [ val(case_info) ]
        ch_target_bed     // channel: [mandatory for WES] [ val(meta), path(bed), path(tbi) ]
        dict          // channel: [ path(dict) ]
        blacklist_bed // channel: [ path(blacklist_bed) ]
        priors        // channel: [ path(priors) ]
        ploidy_model  // channel: [ path(ploidy_model) ]
        cnv_model     // channel: [ path(cnv_model) ]

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

        CALL_CNV_GERMLINECNVCALLER (bam_bai, fasta_no_meta, fai, target_bed, blacklist_bed, dict, priors, ploidy_model, cnv_model )
            .candidate_cnvs_vcf
            .collect{it[1]}
            .set { gcnvcaller_vcf }

        //merge
        tiddit_vcf
            .combine(manta_vcf)
            .toList()
            .set { vcf_list }

        ch_case_info
            .combine(vcf_list)
            .set { merge_input_vcfs }

        SVDB_MERGE (merge_input_vcfs, ["tiddit","manta"])

        TABIX_TABIX (SVDB_MERGE.out.vcf)

        ch_versions = ch_versions.mix(CALL_SV_MANTA.out.versions)
        ch_versions = ch_versions.mix(CALL_SV_TIDDIT.out.versions)
        ch_versions = ch_versions.mix(CALL_CNV_GERMLINECNVCALLER.out.versions)
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:
        vcf      = SVDB_MERGE.out.vcf  // channel: [ val(meta), path(vcf)]
        tbi      = TABIX_TABIX.out.tbi // channel: [ val(meta), path(tbi)]
        versions = ch_versions         // channel: [ path(versions.yml) ]
}
