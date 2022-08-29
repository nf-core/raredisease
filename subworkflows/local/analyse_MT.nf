//
// Analyse MT
//

include { PREPARE_MT_ALIGNMENT } from './prepare_MT_alignment'
include { ALIGN_AND_CALL_MT    } from './align_and_call_MT'



workflow ANALYSE_MT {
    take:
        bam           // channel: [ val(meta), file(bam), file(bai) ]
        index         // channel: [ /path/to/bwamem2/index/ ]
        fasta         // channel: [ genome.fasta ]
        dict          // channel: [ genome.dict ]
        fai           // channel: [ genome.fai ]
        intervals_mt  // channel: [ file(non_control_region.chrM.interval_list) ]

    main:
        ch_versions = Channel.empty()

        // STEP 1: PREPARING MT ALIGNMENT
        PREPARE_MT_ALIGNMENT ( bam )
        ch_versions = ch_versions.mix(PREPARE_MT_ALIGNMENT.out.versions)// Outputs bam files
        
        // STEP 2.1: MT ALLIGNMENT AND VARIANT CALLING
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
    

    emit:
        vcf      = ALIGN_AND_CALL_MT.out.vcf
        tbi      = ALIGN_AND_CALL_MT.out.tbi
        txt      = ALIGN_AND_CALL_MT.out.txt
        html     = ALIGN_AND_CALL_MT.out.html
        versions = ch_versions // channel: [ versions.yml ]
}
