//
// Prepare reference genome files
//

include { BWAMEM2_INDEX } from '../../modules/nf-core/modules/bwamem2/index/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/modules/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY as GATK_SD} from '../../modules/nf-core/modules/gatk4/createsequencedictionary/main'
include { GET_CHROM_SIZES } from '../../modules/local/get_chrom_sizes'

workflow PREPARE_GENOME {
    take:
        fasta // path: genome.fasta

    main:
        ch_fasta    = file(fasta)
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
            ch_fai = SAMTOOLS_FAIDX ( [[], ch_fasta] )
                        .fai
                        .collect{it[1]}
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        }

        ch_sequence_dict = GATK_SD ( ch_fasta ).dict
        ch_chrom_sizes   = GET_CHROM_SIZES ( ch_fai ).sizes
        ch_versions      = ch_versions.mix(GET_CHROM_SIZES.out.versions)


    emit:
        fasta                       = ch_fasta                  // path: genome.fasta
        fai                         = ch_fai                    // path: genome.fasta.fai
        bwamem2_index               = ch_bwamem2_index          // path: bwamem2/index
        chrom_sizes                 = ch_chrom_sizes            // path: chrom.sizes
        sequence_dict               = ch_sequence_dict
        versions                    = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
