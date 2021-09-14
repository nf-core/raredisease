//
// Prepare reference genome files
//


params.bwamem2_idx_options = [:]
params.samtools_faidx_options = [:]

include { BWAMEM2_INDEX } from '../../modules/nf-core/modules/bwamem2/index/main'  addParams( options: params.bwamem2_idx_options )
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/modules/samtools/faidx/main'  addParams( options: params.samtools_faidx_options )

workflow PREPARE_GENOME {
    take:
        fasta // path: genome.fasta

    main:
        ch_fasta = file(fasta)

        ch_bwamem2_version = Channel.empty()
        // Fetch BWAMEM2 index or create from scratch if required
        if ( params.bwamem2 && file(params.bwamem2, checkIfExists:true) ) {
            ch_bwamem2_index = file(params.bwamem2)
        } else {
            ch_bwamem2_index = BWAMEM2_INDEX ( ch_fasta ).index
            ch_bwamem2_version = BWAMEM2_INDEX.out.version
        }

        ch_samtools_version = Channel.empty()
        if ( params.fai ) {
            ch_fai = file(params.fai)
        } else {
            ch_fai = SAMTOOLS_FAIDX ( ch_fasta ).fai
            ch_samtools_version = SAMTOOLS_FAIDX.out.version
        }


    emit:
        fasta                       = ch_fasta                  // path: genome.fasta
        fai                         = ch_fai                    // path: genome.fasta.fai
        samtools_version            = ch_samtools_version       // path: *.version.txt

        bwamem2_index               = ch_bwamem2_index          // path: bwamem2/index
        bwamem2_version             = ch_bwamem2_version        // path: *.version.txt
}
