//
// Prepare reference files
//

include { BWA_INDEX as BWA_INDEX_GENOME                      } from '../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_GENOME              } from '../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_SHIFT_MT            } from '../../modules/nf-core/bwamem2/index/main'
include { CAT_CAT as CAT_CAT_BAIT                            } from '../../modules/nf-core/cat/cat/main'
include { GATK4_BEDTOINTERVALLIST as GATK_BILT               } from '../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_CREATESEQUENCEDICTIONARY as GATK_SD          } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_CREATESEQUENCEDICTIONARY as GATK_SD_SHIFT_MT } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_INTERVALLISTTOOLS as GATK_ILT                } from '../../modules/nf-core/gatk4/intervallisttools/main'
include { GATK4_SHIFTFASTA as GATK_SHIFTFASTA                } from '../../modules/nf-core/gatk4/shiftfasta/main'
include { GET_CHROM_SIZES                                    } from '../../modules/local/get_chrom_sizes'
include { SAMTOOLS_FAIDX as SAMTOOLS_EXTRACT_MT              } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME            } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_MT                } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_MT_SHIFT          } from '../../modules/nf-core/samtools/faidx/main'
include { SENTIEON_BWAINDEX as SENTIEON_BWAINDEX_GENOME      } from '../../modules/local/sentieon/bwamemindex'
include { SENTIEON_BWAINDEX as SENTIEON_BWAINDEX_SHIFT_MT    } from '../../modules/local/sentieon/bwamemindex'
include { TABIX_BGZIPTABIX as TABIX_PBT                      } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_DBSNP                         } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GNOMAD_AF                     } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_PT                            } from '../../modules/nf-core/tabix/tabix/main'
include { UNTAR as UNTAR_VEP_CACHE                           } from '../../modules/nf-core/untar/main'

workflow PREPARE_REFERENCES {
    take:
        ch_fasta           // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai      // channel: [mandatory] [ val(meta), path(fai) ]
        ch_fasta_mt        // channel: [optional for dedicated mt analysis] [ val(meta), path(fasta) ]
        ch_gnomad_af_tab   // channel: [optional; used in for snv annotation] [ val(meta), path(tab) ]
        ch_known_dbsnp     // channel: [optional; used only by sentieon] [ val(meta), path(vcf) ]
        ch_target_bed      // channel: [mandatory for WES] [ path(bed) ]
        ch_vep_cache       // channel: [mandatory for annotation] [ path(cache) ]

    main:
        ch_versions    = Channel.empty()
        ch_tbi         = Channel.empty()
        ch_bgzip_tbi   = Channel.empty()
        ch_bwa         = Channel.empty()
        ch_sentieonbwa = Channel.empty()

        // Genome indices
        BWA_INDEX_GENOME(ch_fasta).index.set{ch_bwa}
        BWAMEM2_INDEX_GENOME(ch_fasta)
        SENTIEON_BWAINDEX_GENOME(ch_fasta).index.set{ch_sentieonbwa}
        SAMTOOLS_FAIDX_GENOME(ch_fasta, [[],[]])
        GATK_SD(ch_fasta)
        GET_CHROM_SIZES( SAMTOOLS_FAIDX_GENOME.out.fai )

        // MT indices
        BWAMEM2_INDEX_SHIFT_MT(ch_fasta_mt)
        SENTIEON_BWAINDEX_SHIFT_MT(ch_fasta_mt)
        ch_fai = Channel.empty().mix(ch_genome_fai, SAMTOOLS_FAIDX_GENOME.out.fai).collect()
        SAMTOOLS_EXTRACT_MT(ch_fasta, ch_fai)
        ch_mt_fasta = Channel.empty().mix(ch_genome_fai, SAMTOOLS_FAIDX_GENOME.out.fai).collect()
        GATK_SD_SHIFT_MT(SAMTOOLS_EXTRACT_MT.out.fa)
        SAMTOOLS_FAIDX_MT_SHIFT(SAMTOOLS_EXTRACT_MT.out.fa, [[],[]])
        GATK_SHIFTFASTA(SAMTOOLS_EXTRACT_MT.out.fa,SAMTOOLS_FAIDX_MT_SHIFT.out.fai, GATK_SD_SHIFT_MT.out.dict)

        // Vcf, tab and bed indices
        TABIX_DBSNP(ch_known_dbsnp)
        TABIX_GNOMAD_AF(ch_gnomad_af_tab)
        TABIX_PT(ch_target_bed).tbi.set { ch_tbi }
        TABIX_PBT(ch_target_bed).gz_tbi.set { ch_bgzip_tbi }

        // Generate bait and target intervals
        GATK_BILT(ch_target_bed, GATK_SD.out.dict).interval_list
        GATK_ILT(GATK_BILT.out.interval_list)
        GATK_ILT.out.interval_list
            .collect{ it[1] }
            .map { it ->
                meta = it[0].toString().split("_split")[0].split("/")[-1] + "_bait.intervals_list"
                return [[id:meta], it]
            }
            .set { ch_bait_intervals_cat_in }
        CAT_CAT_BAIT ( ch_bait_intervals_cat_in )
        UNTAR_VEP_CACHE (ch_vep_cache)

        // Gather versions
        ch_versions = ch_versions.mix(BWA_INDEX_GENOME.out.versions)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX_GENOME.out.versions)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX_SHIFT_MT.out.versions)
        ch_versions = ch_versions.mix(SENTIEON_BWAINDEX_GENOME.out.versions)
        ch_versions = ch_versions.mix(SENTIEON_BWAINDEX_SHIFT_MT.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_GENOME.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_MT_SHIFT.out.versions)
        ch_versions = ch_versions.mix(GATK_SD.out.versions)
        ch_versions = ch_versions.mix(GATK_SD_SHIFT_MT.out.versions)
        ch_versions = ch_versions.mix(GET_CHROM_SIZES.out.versions)
        ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions)
        ch_versions = ch_versions.mix(TABIX_GNOMAD_AF.out.versions)
        ch_versions = ch_versions.mix(TABIX_PT.out.versions)
        ch_versions = ch_versions.mix(TABIX_PBT.out.versions)
        ch_versions = ch_versions.mix(GATK_BILT.out.versions)
        ch_versions = ch_versions.mix(GATK_ILT.out.versions)
        ch_versions = ch_versions.mix(UNTAR_VEP_CACHE.out.versions)

    emit:
        bait_intervals         = CAT_CAT_BAIT.out.file_out.map{ meta, inter -> inter}.collect()   // channel: [ path(intervals) ]
        bwa_index              = Channel.empty().mix(ch_bwa, ch_sentieonbwa).collect()            // channel: [ val(meta), path(index) ]
        bwa_index_mt_shift     = SENTIEON_BWAINDEX_SHIFT_MT.out.index.collect()                   // channel: [ val(meta), path(index) ]
        bwamem2_index          = BWAMEM2_INDEX_GENOME.out.index.collect()                         // channel: [ val(meta), path(index) ]
        bwamem2_index_mt_shift = BWAMEM2_INDEX_SHIFT_MT.out.index.collect()                       // channel: [ val(meta), path(index) ]
        chrom_sizes            = GET_CHROM_SIZES.out.sizes.collect()                              // channel: [ path(sizes) ]
        fai                    = ch_fai                                                           // channel: [ val(meta), path(fai) ]
        fai_mt_shift           = SAMTOOLS_FAIDX_MT_SHIFT.out.fai.collect()                        // channel: [ val(meta), path(fai) ]
        gnomad_af_idx          = TABIX_GNOMAD_AF.out.tbi.collect()                                // channel: [ val(meta), path(fasta) ]
        known_dbsnp_tbi        = TABIX_DBSNP.out.tbi.collect()                                    // channel: [ val(meta), path(fasta) ]
        sequence_dict          = GATK_SD.out.dict.collect()                                       // channel: [ path(dict) ]
        sequence_dict_mt_shift = GATK_SD_SHIFT_MT.out.dict.collect()                              // channel: [ path(dict) ]
        target_bed             = Channel.empty().mix(ch_tbi, ch_bgzip_tbi).collect()              // channel: [ val(meta), path(bed), path(tbi) ]
        target_intervals       = GATK_BILT.out.interval_list.map{ meta, inter -> inter}.collect() // channel: [ path(interval_list) ]
        vep_resources          = UNTAR_VEP_CACHE.out.untar.map{meta, files -> [files]}.collect()  // channel: [ path(cache) ]
        versions               = ch_versions                                                      // channel: [ path(versions.yml) ]

}

