//
// Analyse MT
//
include { PREPARE_MT_ALIGNMENT                           } from './prepare_MT_alignment'
include { ALIGN_AND_CALL_MT                              } from './align_and_call_MT'
include { ALIGN_AND_CALL_MT as ALIGN_AND_CALL_MT_SHIFT   } from './align_and_call_MT'
include { PREPARE_GENOME as PREPARE_GENOME_MT            } from './prepare_genome'
include { PICARD_LIFTOVERVCF                             } from '../../modules/nf-core/picard/liftovervcf/main'

workflow ANALYSE_MT {
    take:
        bam           // channel: [ val(meta), file(bam), file(bai) ]
        index         // channel: [ /path/to/bwamem2/index/ ]
        fasta         // channel: [ genome.fasta ]
        dict          // channel: [ genome.dict ]
        fai           // channel: [ genome.fai ]
        intervals_mt  // channel: [ file(non_control_region.chrM.interval_list) ]
        fasta_shift         // channel: [ genome.fasta ]
        intervals_mt_shift  // channel: [ file(control_region_shifted.chrM.interval_list) ]
        shift_chain

    main:
        ch_versions = Channel.empty()

        // STEP 1: PREPARING MT ALIGNMENT
        PREPARE_MT_ALIGNMENT ( bam )
        ch_versions = ch_versions.mix(PREPARE_MT_ALIGNMENT.out.versions)// Outputs bam files

        //STEP 2.1: MT ALLIGNMENT  AND VARIANT CALLING
        ch_intervals_mt = Channel.fromPath(params.intervals_mt)
        ALIGN_AND_CALL_MT (
            PREPARE_MT_ALIGNMENT.out.fastq,
            PREPARE_MT_ALIGNMENT.out.bam,
            index,
            fasta,
            dict,
            fai,
            ch_intervals_mt
            )
        ch_versions = ch_versions.mix(ALIGN_AND_CALL_MT.out.versions)

        // STEP 2.2: MT ALLIGNMENT SHIFT AND VARIANT CALLING
        ch_intervals_mt_shift = Channel.fromPath(params.intervals_mt_shift)
        PREPARE_GENOME_MT("bwamem2",[],[],fasta_shift ,[],[],[],false).set { ch_genome }
        ch_versions = ch_versions.mix(ch_genome.versions)
        ch_dict_shift = ch_genome.sequence_dict
        ch_fai_shift = ch_genome.fai
        ch_index_shift =ch_genome.bwamem2_index

        ALIGN_AND_CALL_MT_SHIFT (
            PREPARE_MT_ALIGNMENT.out.fastq,
            PREPARE_MT_ALIGNMENT.out.bam,
            ch_index_shift,
            fasta_shift,
            ch_dict_shift,
            ch_fai_shift,
            ch_intervals_mt_shift
            )
        ch_versions = ch_versions.mix(ALIGN_AND_CALL_MT_SHIFT.out.versions)

        // STEP 2.3: PICARD_LIFTOVERVCF
        PICARD_LIFTOVERVCF (
              ALIGN_AND_CALL_MT_SHIFT.out.vcf,
              dict,
              shift_chain,
              fasta)
        ch_versions = ch_versions.mix(PICARD_LIFTOVERVCF.out.versions)


    emit:
        vcf         = ALIGN_AND_CALL_MT.out.vcf
        tbi         = ALIGN_AND_CALL_MT.out.tbi
        txt         = ALIGN_AND_CALL_MT.out.txt
        html        = ALIGN_AND_CALL_MT.out.html
        vcf_shift   = ALIGN_AND_CALL_MT_SHIFT.out.vcf
        tbi_shift   = ALIGN_AND_CALL_MT_SHIFT.out.tbi
        txt_shift   = ALIGN_AND_CALL_MT_SHIFT.out.txt
        html_shift  = ALIGN_AND_CALL_MT_SHIFT.out.html
        vcf_lift    = PICARD_LIFTOVERVCF.out.vcf_lifted
        vcf_unlift  = PICARD_LIFTOVERVCF.out.vcf_unlifted
        versions    = ch_versions // channel: [ versions.yml ]

}
