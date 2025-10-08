//
// Prepare reference files
//

include { BWA_INDEX as BWA_INDEX_GENOME                      } from '../../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_GENOME              } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEME_INDEX as BWAMEME_INDEX_GENOME              } from '../../../modules/nf-core/bwameme/index/main'
include { CAT_CAT as CAT_CAT_BAIT                            } from '../../../modules/nf-core/cat/cat/main'
include { GATK4_CREATESEQUENCEDICTIONARY as GATK_SD          } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GET_CHROM_SIZES                                    } from '../../../modules/local/get_chrom_sizes'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME            } from '../../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_REFERENCES {
    take:
        ch_genome_fasta              // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai                // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_dictionary         // channel: [mandatory] [ val(meta), path(fai) ]

    main:
        ch_versions      = Channel.empty()
        ch_tbi           = Channel.empty()
        ch_bgzip_tbi     = Channel.empty()
        ch_bwa           = Channel.empty()
        ch_sentieonbwa   = Channel.empty()

        // Genome indices
        SAMTOOLS_FAIDX_GENOME(ch_genome_fasta, [[],[]])
        GATK_SD(ch_genome_fasta)
        ch_fai  = Channel.empty().mix(ch_genome_fai, SAMTOOLS_FAIDX_GENOME.out.fai).collect()
        ch_dict = Channel.empty().mix(ch_genome_dictionary, GATK_SD.out.dict).collect()
        GET_CHROM_SIZES( ch_fai )

        // Genome alignment indices
        BWA_INDEX_GENOME(ch_genome_fasta).index.set{ch_bwa}
        BWAMEM2_INDEX_GENOME(ch_genome_fasta)
        BWAMEME_INDEX_GENOME(ch_genome_fasta)

        // Gather versions
        ch_versions = ch_versions.mix(BWA_INDEX_GENOME.out.versions)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX_GENOME.out.versions)
        ch_versions = ch_versions.mix(BWAMEME_INDEX_GENOME.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_GENOME.out.versions)

    emit:
        genome_bwa_index      = Channel.empty().mix(ch_bwa, ch_sentieonbwa).collect()                        // channel: [ val(meta), path(index) ]
        genome_bwamem2_index  = BWAMEM2_INDEX_GENOME.out.index.collect()                                     // channel: [ val(meta), path(index) ]
        genome_bwameme_index  = BWAMEME_INDEX_GENOME.out.index.collect()                                     // channel: [ val(meta), path(index) ]
        genome_chrom_sizes    = GET_CHROM_SIZES.out.sizes.collect()                                          // channel: [ path(sizes) ]
        genome_fai            = ch_fai                                                                       // channel: [ val(meta), path(fai) ]
        genome_dict           = ch_dict                                                                      // channel: [ val(meta), path(dict) ]
        versions              = ch_versions                                                                  // channel: [ path(versions.yml) ]

}
