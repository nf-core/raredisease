//
// Prepare reference genome files
//


// params.bwamem2_idx_options = [:]
// params.samtools_faidx_options = [:]

include { BWAMEM2_INDEX } from '../../modules/nf-core/modules/bwamem2/index/main'  //addParams( options: params.bwamem2_idx_options )
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/modules/samtools/faidx/main'  //addParams( options: params.samtools_faidx_options )

workflow PREPARE_GENOME {
    take:
        fasta // path: genome.fasta

    main:
        ch_fasta = file(fasta)
        ch_versions = Channel.empty()

        // Fetch BWAMEM2 index or create from scratch if required
        if ( params.aligner == 'bwamem2' ) {
            if ( params.bwamem2_index && file(params.bwamem2_index, checkIfExists:true) ) {
                ch_bwamem2_index = file(params.bwamem2_index)
            } else {
                ch_bwamem2_index = BWAMEM2_INDEX ( ch_fasta ).index
                ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
            }
        }

        if ( params.fasta_fai ) {
            ch_fai = file(params.fasta_fai)
        } else {
            ch_fai = SAMTOOLS_FAIDX ( ch_fasta ).fai
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        }


    emit:
        fasta                       = ch_fasta                  // path: genome.fasta
        fai                         = ch_fai                    // path: genome.fasta.fai
        bwamem2_index               = ch_bwamem2_index          // path: bwamem2/index

        versions                    = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
