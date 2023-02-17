//
// Prepare reference files
//

include { BWA_INDEX as BWA_INDEX_GENOME                      } from '../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_GENOME              } from '../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_SHIFT_MT            } from '../../modules/nf-core/bwamem2/index/main'
include { CAT_CAT as CAT_CAT_BAIT                            } from '../../modules/nf-core/cat/cat/main'
include { CHECK_VCF                                          } from './preprocessing/prepare_vcf'
include { GATK4_BEDTOINTERVALLIST as GATK_BILT               } from '../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_CREATESEQUENCEDICTIONARY as GATK_SD          } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_CREATESEQUENCEDICTIONARY as GATK_SD_SHIFT_MT } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_INTERVALLISTTOOLS as GATK_ILT                } from '../../modules/nf-core/gatk4/intervallisttools/main'
include { GET_CHROM_SIZES                                    } from '../../modules/local/get_chrom_sizes'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME            } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SHIFT_MT          } from '../../modules/nf-core/samtools/faidx/main'
include { SENTIEON_BWAINDEX as SENTIEON_BWAINDEX_GENOME      } from '../../modules/local/sentieon/bwamemindex'
include { SENTIEON_BWAINDEX as SENTIEON_BWAINDEX_SHIFT_MT    } from '../../modules/local/sentieon/bwamemindex'
include { TABIX_BGZIPTABIX as TABIX_PBT                      } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_DBSNP                         } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GNOMAD_AF                     } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_PT                            } from '../../modules/nf-core/tabix/tabix/main'
include { UNTAR as UNTAR_VEP_CACHE                           } from '../../modules/nf-core/untar/main'

workflow PREPARE_REFERENCES {
    take:
        fasta_no_meta               // [mandatory] genome.fasta
        fasta_meta
        fai                 // [optional ] genome.fai
        mt_fasta_shift_no_meta
        mt_fasta_shift_meta
        gnomad_af_tab
        gnomad_vcf_in
        known_dbsnp
        target_bed
        vep_cache

    main:
        ch_versions    = Channel.empty()
        ch_tbi         = Channel.empty()
        ch_bgzip_tbi   = Channel.empty()
        ch_bwa         = Channel.empty()
        ch_sentieonbwa = Channel.empty()

        // Genome indices
        BWA_INDEX_GENOME(fasta_meta).index.set{ch_bwa}
        BWAMEM2_INDEX_GENOME(fasta_meta)
        BWAMEM2_INDEX_SHIFT_MT(mt_fasta_shift_meta)
        SENTIEON_BWAINDEX_GENOME(fasta_meta).index.set{ch_sentieonbwa}
        SENTIEON_BWAINDEX_SHIFT_MT(mt_fasta_shift_meta)
        SAMTOOLS_FAIDX_GENOME(fasta_meta)
        SAMTOOLS_FAIDX_SHIFT_MT(mt_fasta_shift_meta)
        GATK_SD(fasta_no_meta)
        GATK_SD_SHIFT_MT(mt_fasta_shift_no_meta)
        GET_CHROM_SIZES( SAMTOOLS_FAIDX_GENOME.out.fai )

        // Vcf, tab and bed indices
        TABIX_DBSNP(known_dbsnp)
        TABIX_GNOMAD_AF(gnomad_af_tab)
        TABIX_PT(target_bed).tbi.set { ch_tbi }
        TABIX_PBT(target_bed).gz_tbi.set { ch_bgzip_tbi }

        // Check if a vcf file is normalized and create index
        CHECK_VCF(gnomad_vcf_in, fasta_no_meta)

        // Generate bait and target intervals
        GATK_BILT(target_bed, GATK_SD.out.dict).interval_list
        GATK_ILT(GATK_BILT.out.interval_list)
        GATK_ILT.out.interval_list
            .collect{ it[1] }
            .map { it ->
                meta = it[0].toString().split("_split")[0].split("/")[-1] + "_bait.intervals_list"
                return [[id:meta], it]
            }
            .set { ch_bait_intervals_cat_in }
        CAT_CAT_BAIT ( ch_bait_intervals_cat_in )
        UNTAR_VEP_CACHE (vep_cache)

        // Gather versions
        ch_versions = ch_versions.mix(BWA_INDEX_GENOME.out.versions)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX_GENOME.out.versions)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX_SHIFT_MT.out.versions)
        ch_versions = ch_versions.mix(SENTIEON_BWAINDEX_GENOME.out.versions)
        ch_versions = ch_versions.mix(SENTIEON_BWAINDEX_SHIFT_MT.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_GENOME.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_SHIFT_MT.out.versions)
        ch_versions = ch_versions.mix(GATK_SD.out.versions)
        ch_versions = ch_versions.mix(GATK_SD_SHIFT_MT.out.versions)
        ch_versions = ch_versions.mix(GET_CHROM_SIZES.out.versions)
        ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions)
        ch_versions = ch_versions.mix(TABIX_GNOMAD_AF.out.versions)
        ch_versions = ch_versions.mix(CHECK_VCF.out.versions)
        ch_versions = ch_versions.mix(TABIX_PT.out.versions)
        ch_versions = ch_versions.mix(TABIX_PBT.out.versions)
        ch_versions = ch_versions.mix(GATK_BILT.out.versions)
        ch_versions = ch_versions.mix(GATK_ILT.out.versions)
        ch_versions = ch_versions.mix(UNTAR_VEP_CACHE.out.versions)

    emit:
        bait_intervals            = CAT_CAT_BAIT.out.file_out.map { id, it -> [it] }.collect()
        bwa_index                 = Channel.empty().mix(ch_bwa, ch_sentieonbwa).collect()
        bwa_index_mt_shift        = SENTIEON_BWAINDEX_SHIFT_MT.out.index.collect()
        bwamem2_index             = BWAMEM2_INDEX_GENOME.out.index.collect()
        bwamem2_index_mt_shift    = BWAMEM2_INDEX_SHIFT_MT.out.index.collect()
        chrom_sizes               = GET_CHROM_SIZES.out.sizes.collect()
        fasta_fai                 = SAMTOOLS_FAIDX_GENOME.out.fai.map{ meta, fai -> [fai] }.collect()
        fasta_fai_meta            = SAMTOOLS_FAIDX_GENOME.out.fai.collect()
        fasta_fai_mt_shift        = SAMTOOLS_FAIDX_SHIFT_MT.out.fai.map{ meta, fai -> [fai] }.collect()
        gnomad_af_idx             = TABIX_GNOMAD_AF.out.tbi.collect()
        gnomad_tbi                = CHECK_VCF.out.index.collect()
        gnomad_vcf                = CHECK_VCF.out.vcf.collect()
        known_dbsnp_tbi           = TABIX_DBSNP.out.tbi.collect()
        sequence_dict             = GATK_SD.out.dict.collect()
        sequence_dict_meta        = GATK_SD.out.dict.map {it -> [[id:it[0].simpleName], it]}.collect()
        sequence_dict_mt_shift    = GATK_SD_SHIFT_MT.out.dict.collect()
        target_bed                = Channel.empty().mix(ch_tbi, ch_bgzip_tbi).collect()
        target_intervals          = GATK_BILT.out.interval_list.collect{it[1]}.collect()
        vep_resources             = UNTAR_VEP_CACHE.out.untar.map{meta, files -> [files]}collect()
        versions                  = ch_versions

}

