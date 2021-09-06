//
// Prepare indices for reference genome files
//

params.bwamem2_idx_options = [:]

include { BWAMEM2_INDEX } from '../../modules/nf-core/modules/bwamem2/index/main'  addParams( options: params.bwamem2_idx_options )

workflow PREPARE_GENOME {
    take:
        prepare_tool_indicies // list: tools to prepare indices for

    main:

        ch_bwamem2_index = Channel.empty()
        ch_bwamem2_version = Channel.empty()
        // Fetch BWAMEM2 index or create from scratch if required
        if ('bwamem2' in prepare_tool_indicies) {
            if (params.bwamem2_index) {
                ch_bwamem2_index = file(params.bwamem2_index)
            } else if (params.bwamem2){
                ch_bwamem2_index = file(params.bwamem2)
            } else {
                ch_bwamem2_index = BWAMEM2_INDEX ( params.fasta ).index
                ch_bwamem2_version = BWAMEM2_INDEX.out.version
            }
        }


    emit:
        bwamem2_index = ch_bwamem2_index
        bwamem2_version = ch_bwamem2_version
}
