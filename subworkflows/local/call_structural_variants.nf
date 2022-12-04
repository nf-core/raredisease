//
// A nested subworkflow to call structural variants.
//

include { CALL_SV_MANTA     } from './variant_calling/call_sv_manta'
include { CALL_SV_TIDDIT    } from './variant_calling/call_sv_tiddit'
include { SVDB_MERGE        } from '../../modules/nf-core/svdb/merge/main'
include { CALL_CNV_CNVPYTOR } from './variant_calling/call_cnv_cnvpytor'

workflow CALL_STRUCTURAL_VARIANTS {

    take:
        bam           // channel: [ val(meta), path(bam) ]
        bai           // channel: [ val(meta), path(bai) ]
        bam_bai       // channel: [ val(meta), path(bam), path(bai) ]
        bwa_index     // channel: [ val(meta), path(bwa_index)]
        fasta_no_meta // channel: [ path(genome.fasta) ]
        fasta_meta    // channel: [ val(meta), path(genome.fasta) ]
        fai           // channel: [ path(genome.fai) ]
        case_info     // channel: [ val(case_info) ]
        target_bed    // channel: [ path(target.bed) ]
        cnvpytor_bins // channel: [ val("binsizes") ]

    main:
        ch_versions = Channel.empty()

        CALL_SV_MANTA (bam, bai, fasta_no_meta, fai, case_info, target_bed)
            .diploid_sv_vcf
            .collect{it[1]}
            .set{ manta_vcf }

        CALL_SV_TIDDIT (bam_bai, fasta_meta, bwa_index, case_info)
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

        case_info
            .combine(vcf_list)
            .set { merge_input_vcfs }

        // SVDB_MERGE ( merge_input_vcfs, ["tiddit","manta","cnvpytor"] )
        SVDB_MERGE (merge_input_vcfs, ["tiddit","manta"])

        ch_versions = ch_versions.mix(CALL_SV_MANTA.out.versions)
        ch_versions = ch_versions.mix(CALL_SV_TIDDIT.out.versions)
        // ch_versions = ch_versions.mix(CALL_CNV_CNVPYTOR.out.versions)

    emit:
        vcf        = SVDB_MERGE.out.vcf
        versions   = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}
