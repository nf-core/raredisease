//
// Analyse MT
//
include { CONVERT_MT_BAM_TO_FASTQ                        } from './mitochondria/convert_mt_bam_to_fastq'
include { ALIGN_AND_CALL_MT                              } from './mitochondria/align_and_call_MT'
include { ALIGN_AND_CALL_MT as ALIGN_AND_CALL_MT_SHIFT   } from './mitochondria/align_and_call_MT'
include { PICARD_LIFTOVERVCF                             } from '../../modules/nf-core/picard/liftovervcf/main'
include { MERGE_ANNOTATE_MT                              } from './mitochondria/merge_annotate_MT'

workflow ANALYSE_MT {
    take:
        bam                    // channel: [ val(meta), file(bam), file(bai) ]
        genome_bwamem2_index   // channel: [ /path/to/bwamem2/index/ ]
        genome_fasta           // channel: [ genome.fasta ]
        genome_dict            // channel: [ genome.dict ]
        genome_fai             // channel: [ genome.fai ]
        mt_intervals           // channel: [ file(non_control_region.chrM.interval_list) ]
        shift_mt_bwamem2_index // channel: [ /path/to/bwamem2/index/ ]
        shift_mt_fasta         // channel: [ genome.fasta ]
        shift_mt_dict          // channel: [ genome.dict ]
        shift_mt_fai           // channel: [ genome.fai ]
        shift_mt_intervals     // channel: [ file(control_region_shifted.chrM.interval_list) ]
        shift_mt_backchain     // channel: [ file(shift.back_chain) ]
        vep_genome
        vep_cache_version
        vep_cache
        case_info              // channel: [ val(case_info) ]

    main:
        ch_versions = Channel.empty()

        // STEP 1: PREPARING MT ALIGNMENT
        CONVERT_MT_BAM_TO_FASTQ(bam)
        ch_versions = ch_versions.mix(CONVERT_MT_BAM_TO_FASTQ.out.versions)// Outputs bam files

        //STEP 2.1: MT ALLIGNMENT  AND VARIANT CALLING
        ALIGN_AND_CALL_MT (
            CONVERT_MT_BAM_TO_FASTQ.out.fastq,
            CONVERT_MT_BAM_TO_FASTQ.out.bam,
            genome_bwamem2_index,
            genome_fasta,
            genome_dict,
            genome_fai,
            mt_intervals
            )
        ch_versions = ch_versions.mix(ALIGN_AND_CALL_MT.out.versions)

        ALIGN_AND_CALL_MT_SHIFT (
            CONVERT_MT_BAM_TO_FASTQ.out.fastq,
            CONVERT_MT_BAM_TO_FASTQ.out.bam,
            shift_mt_bwamem2_index,
            shift_mt_fasta,
            shift_mt_dict,
            shift_mt_fai,
            shift_mt_intervals
            )
        ch_versions = ch_versions.mix(ALIGN_AND_CALL_MT_SHIFT.out.versions)

        // STEP 2.3: PICARD_LIFTOVERVCF
        PICARD_LIFTOVERVCF (
            ALIGN_AND_CALL_MT_SHIFT.out.vcf,
            genome_dict,
            shift_mt_backchain,
            genome_fasta)
        ch_versions = ch_versions.mix(PICARD_LIFTOVERVCF.out.versions)

        // STEP 3: MT MERGE AND ANNOTATE VARIANTS
        MERGE_ANNOTATE_MT(
            ALIGN_AND_CALL_MT.out.vcf,
            PICARD_LIFTOVERVCF.out.vcf_lifted,
            genome_fasta,
            genome_dict,
            genome_fai,
            vep_genome,
            vep_cache_version,
            vep_cache,
            case_info)
        ch_versions = ch_versions.mix(MERGE_ANNOTATE_MT.out.versions)

    emit:
        vcf          = MERGE_ANNOTATE_MT.out.vcf
        tbi          = MERGE_ANNOTATE_MT.out.tbi
        stats        = ALIGN_AND_CALL_MT.out.stats
        filt_sats    = ALIGN_AND_CALL_MT.out.filt_sats
        stats_sh     = ALIGN_AND_CALL_MT_SHIFT.out.stats
        filt_sats_sh = ALIGN_AND_CALL_MT_SHIFT.out.filt_sats
        haplog       = MERGE_ANNOTATE_MT.out.haplog
        report       = MERGE_ANNOTATE_MT.out.report
        txt          = ALIGN_AND_CALL_MT.out.txt
        html         = ALIGN_AND_CALL_MT.out.html
        txt_sh       = ALIGN_AND_CALL_MT_SHIFT.out.txt
        html_sh      = ALIGN_AND_CALL_MT_SHIFT.out.html
        versions     = ch_versions // channel: [ versions.yml ]
}
