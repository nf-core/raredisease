#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BWAMEM2_INDEX        } from '../../../../modules/nf-core/modules/bwamem2/index/main'

include { PREPARE_MT_ALIGNMENT } from '../../../../subworkflows/local/prepare_MT_alignment'
include { ALIGN_MT             } from '../../../../subworkflows/local/align_MT'




workflow align_mt {
    input        = [ [ id:'A', sample:'NA12878', single_end: false], // meta map
                    file(params.test_data['mini_human_genome']['alignment']['bam'], checkIfExists: true),
                    file(params.test_data['mini_human_genome']['alignment']['bai'], checkIfExists: true),
                    ]

    fasta        = [ file(params.test_data['mini_human_genome']['fasta'], checkIfExists: true) ]
    fai          = [ file(params.test_data['mini_human_genome']['fai'], checkIfExists: true) ]
    dict         = [ file(params.test_data['mini_human_genome']['dict'], checkIfExists: true) ]
    intervals_mt = [ file(params.test_data['mitochondria']['intervals_mt'], checkIfExists: true ) ]

    PREPARE_MT_ALIGNMENT ( input )

    BWAMEM2_INDEX ( fasta )

    ALIGN_MT (
        PREPARE_MT_ALIGNMENT.out.fastq,
        PREPARE_MT_ALIGNMENT.out.bam,
        BWAMEM2_INDEX.out.index,
        fasta,
        dict,
        fai,
        intervals_mt
    )
}
