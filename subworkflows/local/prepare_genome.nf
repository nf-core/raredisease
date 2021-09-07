//
// Prepare indices for reference genome files
//

params.bwamem2_idx_options = [:]

include { BWAMEM2_INDEX } from '../../modules/nf-core/modules/bwamem2/index/main'  addParams( options: params.bwamem2_idx_options )

workflow PREPARE_GENOME {
    take:
        fasta

    main:

        ch_bwamem2_index = Channel.empty()
        ch_bwamem2_version = Channel.empty()
        // Fetch BWAMEM2 index or create from scratch if required
        if ( params.bwamem2 && file(params.bwamem2, checkIfExists:true) ) {
            ch_bwamem2_index = file(params.bwamem2)
        } else {
            ch_bwamem2_index = BWAMEM2_INDEX ( fasta ).index
            ch_bwamem2_version = BWAMEM2_INDEX.out.version
        }


    emit:
        bwamem2_index = ch_bwamem2_index
        bwamem2_version = ch_bwamem2_version
}
