//
// Prepare reference files
//

include { BEDTOOLS_SLOP as BEDTOOLS_PAD_TARGET_BED           } from '../../../modules/nf-core/bedtools/slop/main'
include { BWA_INDEX as BWA_INDEX_GENOME                      } from '../../../modules/nf-core/bwa/index/main'
include { BWA_INDEX as BWA_INDEX_MT                          } from '../../../modules/nf-core/bwa/index/main'
include { BWA_INDEX as BWA_INDEX_MT_SHIFT                    } from '../../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_GENOME              } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_MT                  } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_MT_SHIFT            } from '../../../modules/nf-core/bwamem2/index/main'
include { BWAMEME_INDEX as BWAMEME_INDEX_GENOME              } from '../../../modules/nf-core/bwameme/index/main'
include { CAT_CAT as CAT_CAT_BAIT                            } from '../../../modules/nf-core/cat/cat/main'
include { GATK4_BEDTOINTERVALLIST as GATK_BILT               } from '../../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_CREATESEQUENCEDICTIONARY as GATK_SD          } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_CREATESEQUENCEDICTIONARY as GATK_SD_MT       } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_INTERVALLISTTOOLS as GATK_ILT                } from '../../../modules/nf-core/gatk4/intervallisttools/main'
include { GATK4_SHIFTFASTA as GATK_SHIFTFASTA                } from '../../../modules/nf-core/gatk4/shiftfasta/main'
include { GET_CHROM_SIZES                                    } from '../../../modules/local/get_chrom_sizes'
include { RTGTOOLS_FORMAT                                    } from '../../../modules/nf-core/rtgtools/format/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_EXTRACT_MT              } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME            } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_MT                } from '../../../modules/nf-core/samtools/faidx/main'
include { SENTIEON_BWAINDEX as SENTIEON_BWAINDEX_GENOME      } from '../../../modules/nf-core/sentieon/bwaindex/main'
include { SENTIEON_BWAINDEX as SENTIEON_BWAINDEX_MT          } from '../../../modules/nf-core/sentieon/bwaindex/main'
include { SENTIEON_BWAINDEX as SENTIEON_BWAINDEX_MT_SHIFT    } from '../../../modules/nf-core/sentieon/bwaindex/main'
include { TABIX_BGZIPTABIX as TABIX_PBT                      } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPINDEX_PADDED_BED    } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPINDEX_VCFANNOEXTRA  } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_VCFANNOEXTRA                  } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_DBSNP                         } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GNOMAD_AF                     } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_PT                            } from '../../../modules/nf-core/tabix/tabix/main'
include { UNTAR as UNTAR_VEP_CACHE                           } from '../../../modules/nf-core/untar/main'

workflow PREPARE_REFERENCES {
    take:
        ch_genome_fasta              // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_gnomad_af_tab             // channel: [optional; used in for snv annotation] [ val(meta), path(tab) ]
        ch_known_dbsnp               // channel: [optional; used only by sentieon] [ val(meta), path(vcf) ]
        ch_target_bed                // channel: [mandatory for WES] [ path(bed) ]
        ch_vcfanno_extra_unprocessed // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_vep_cache                 // channel: [mandatory for annotation] [ path(cache) ]
        val_aligner                  // String: "bwa", "bwamem2", "sentieon" or "bwameme"
        val_analysis_type            // String: "wgs", "wes", or "mito"
        val_bwa                      // String: path to bwa index
        val_bwamem2                  // String: path to bwamem2 index
        val_bwameme                  // String: path to bwameme index
        val_fai                      // String: path to genome fasta index
        val_mtaligner                // String: "bwa", "bwamem2", or "sentieon"
        val_mtfasta                  // String: path to mitochondrial fasta
        val_run_mt_for_wes           // Boolean
        val_genome_dict              // String: path to genome dictionary

    main:
        ch_versions              = channel.empty()
        ch_bgzip_tbi             = channel.empty()
        ch_bwa                   = channel.empty()
        ch_genome_bwameme_index  = channel.empty()
        ch_genome_bwamem2_index  = channel.empty()
        ch_mt_bwa_index          = channel.empty()
        ch_mt_bwamem2_index      = channel.empty()
        ch_mtshift_bwa_index     = channel.empty()
        ch_mtshift_bwamem2_index = channel.empty()
        ch_mt_dict               = channel.empty()
        ch_mt_fai                = channel.empty()
        ch_mt_fasta              = channel.empty()
        ch_sentieonbwa           = channel.empty()
        ch_vcfanno_extra         = channel.empty()
        ch_vcfanno_bgzip         = channel.empty()
        ch_vcfanno_index         = channel.empty()

        // Genome indices
        if (!val_fai) {
            ch_fai = SAMTOOLS_FAIDX_GENOME(ch_genome_fasta, [[],[]]).fai.collect()
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_GENOME.out.versions)
        } else {
            ch_fai  = channel.fromPath(val_fai).map {it -> [[id:it.simpleName], it]}.collect()
        }

        if (!val_genome_dict) {
            ch_genome_dict = GATK_SD(ch_genome_fasta).dict.collect()
            ch_versions = ch_versions.mix(GATK_SD.out.versions)
        } else {
            ch_genome_dict = channel.fromPath(val_genome_dict).map {it -> [[id:it.simpleName], it]}.collect()
        }

        GET_CHROM_SIZES( ch_fai )

        // Genome alignment indices
        if (!val_bwa) {
            if (!val_aligner.equals("sentieon") || val_mtaligner.equals("bwa")) {
                BWA_INDEX_GENOME(ch_genome_fasta).index.set{ch_bwa}
                ch_versions = ch_versions.mix(BWA_INDEX_GENOME.out.versions)
            }
            if (val_aligner.equals("sentieon") || val_mtaligner.equals("sentieon")) {
                SENTIEON_BWAINDEX_GENOME(ch_genome_fasta).index.set{ch_sentieonbwa}
                ch_versions = ch_versions.mix(SENTIEON_BWAINDEX_GENOME.out.versions)
            }
        }
        if (!val_bwamem2 && (val_aligner.equals("bwamem2") || val_mtaligner.equals("bwamem2"))) {
            ch_genome_bwamem2_index = BWAMEM2_INDEX_GENOME(ch_genome_fasta).index.collect()
            ch_versions = ch_versions.mix(BWAMEM2_INDEX_GENOME.out.versions)
        }
        if (!val_bwamem2 && val_aligner.equals("bwameme")) {
            ch_genome_bwameme_index = BWAMEME_INDEX_GENOME(ch_genome_fasta).index.collect()
            ch_versions = ch_versions.mix(BWAMEME_INDEX_GENOME.out.versions)
        }

        // MT genome indices

        if (val_analysis_type.matches("wgs|mito") || val_run_mt_for_wes) {
            if (!val_mtfasta) {
                ch_mt_fasta = SAMTOOLS_EXTRACT_MT(ch_genome_fasta, ch_fai).fa.collect()
                ch_versions = ch_versions.mix(SAMTOOLS_EXTRACT_MT.out.versions)
            } else {
                ch_mt_fasta = channel.fromPath(val_mtfasta).map { it -> [[id:it.simpleName], it] }.collect()
            }

            ch_mt_fai  = SAMTOOLS_FAIDX_MT(ch_mt_fasta, [[],[]]).fai.collect()
            ch_mt_dict = GATK_SD_MT(ch_mt_fasta).dict.collect()
            GATK_SHIFTFASTA(ch_mt_fasta, ch_mt_fai, ch_mt_dict)

            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_MT.out.versions,GATK_SD_MT.out.versions, GATK_SHIFTFASTA.out.versions)
        }

        // MT alignment indices
        if ((val_analysis_type.matches("wgs|mito") || val_run_mt_for_wes) && val_mtaligner.equals("bwamem2")) {
            ch_mt_bwamem2_index      = BWAMEM2_INDEX_MT(ch_mt_fasta).index.collect()
            ch_mtshift_bwamem2_index = BWAMEM2_INDEX_MT_SHIFT(GATK_SHIFTFASTA.out.shift_fa).index.collect()
            ch_versions              = ch_versions.mix(BWAMEM2_INDEX_MT.out.versions, BWAMEM2_INDEX_MT_SHIFT.out.versions)
        }
        if ((val_analysis_type.matches("wgs|mito") || val_run_mt_for_wes) && val_mtaligner.equals("bwa")) {
            ch_mt_bwa_index          = BWA_INDEX_MT(ch_mt_fasta).index.collect()
            ch_mtshift_bwa_index     = BWA_INDEX_MT_SHIFT(GATK_SHIFTFASTA.out.shift_fa).index.collect()
            ch_versions              = ch_versions.mix(BWA_INDEX_MT.out.versions, BWA_INDEX_MT_SHIFT.out.versions)
        }
        if ((val_analysis_type.matches("wgs|mito") || val_run_mt_for_wes) && val_mtaligner.equals("sentieon")) {
            ch_mt_bwa_index          = SENTIEON_BWAINDEX_MT(ch_mt_fasta).index.collect()
            ch_mtshift_bwa_index     = SENTIEON_BWAINDEX_MT_SHIFT(GATK_SHIFTFASTA.out.shift_fa).index.collect()
            ch_versions              = ch_versions.mix(SENTIEON_BWAINDEX_MT.out.versions, SENTIEON_BWAINDEX_MT_SHIFT.out.versions)
        }

        GATK_SHIFTFASTA.out.intervals
            .multiMap{ _meta, files ->
                    shift_intervals:
                        def ind = files.findIndexValues {file -> file.toString().endsWith("shifted.intervals")}
                        files[ind]
                    intervals:
                        ind = files.findIndexValues {file -> !(file.toString().endsWith("shifted.intervals"))}
                        files[ind]
            }
            .set {ch_shiftfasta_mtintervals}

        // Vcf, tab and bed indices
        TABIX_DBSNP(ch_known_dbsnp)
        TABIX_GNOMAD_AF(ch_gnomad_af_tab)

        // Index target bed file in case of gz input
        TABIX_PT(ch_target_bed)
        ch_target_bed
            .join(TABIX_PT.out.tbi)
            .set{ ch_trgt_bed_tbi }
        // Compress and index target bed file in case of uncompressed input
        TABIX_PBT(ch_target_bed).gz_tbi
            .set { ch_bgzip_tbi }
        ch_target_bed_gz_tbi = channel.empty()
            .mix(ch_trgt_bed_tbi, ch_bgzip_tbi)

        ch_vcfanno_extra_unprocessed
            .branch { _meta, vcf ->
                bgzipindex: !vcf.toString().endsWith(".gz")
                index: vcf.toString().endsWith(".gz")
            }
            .set { ch_vcfanno_tabix_in }

        TABIX_VCFANNOEXTRA(ch_vcfanno_tabix_in.index).tbi
            .join(ch_vcfanno_tabix_in.index)
            .map { _meta, tbi, vcf -> return [[vcf,tbi]]}
            .set {ch_vcfanno_index}

        TABIX_BGZIPINDEX_VCFANNOEXTRA(ch_vcfanno_tabix_in.bgzipindex)
        channel.empty()
            .mix(TABIX_BGZIPINDEX_VCFANNOEXTRA.out.gz_tbi, TABIX_BGZIPINDEX_VCFANNOEXTRA.out.gz_csi)
            .map { _meta, vcf, index -> return [[vcf,index]] }
            .set {ch_vcfanno_bgzip}

        channel.empty()
            .mix(ch_vcfanno_bgzip, ch_vcfanno_index)
            .collect()
            .set{ch_vcfanno_extra}

        // Pad bed file
        BEDTOOLS_PAD_TARGET_BED(
            ch_target_bed,
            ch_fai.map { _meta, fai -> return fai }
        )
        TABIX_BGZIPINDEX_PADDED_BED(BEDTOOLS_PAD_TARGET_BED.out.bed).gz_tbi
            .set { ch_target_bed_gz_tbi }

        // Generate bait and target intervals
        GATK_BILT(ch_target_bed, ch_genome_dict).interval_list
        GATK_ILT(GATK_BILT.out.interval_list)
        GATK_ILT.out.interval_list
            .collect{ _meta, list -> list }
            .map { list ->
                def meta = list.toString().split("_split")[0].split("/")[-1] + "_bait.intervals_list"
                return [[id:meta], list]
            }
            .set { ch_bait_intervals_cat_in }
        CAT_CAT_BAIT ( ch_bait_intervals_cat_in )
        UNTAR_VEP_CACHE (ch_vep_cache)

        // RTG tools
        ch_genome_fasta.map { meta, fasta -> return [meta, fasta, [], [] ] }
            .set {ch_rtgformat_in}
        RTGTOOLS_FORMAT(ch_rtgformat_in)

        // Gather versions
        ch_versions = ch_versions.mix(GET_CHROM_SIZES.out.versions)
        ch_versions = ch_versions.mix(TABIX_GNOMAD_AF.out.versions)
        ch_versions = ch_versions.mix(TABIX_PT.out.versions)
        ch_versions = ch_versions.mix(TABIX_PBT.out.versions)
        ch_versions = ch_versions.mix(TABIX_BGZIPINDEX_VCFANNOEXTRA.out.versions)
        ch_versions = ch_versions.mix(TABIX_VCFANNOEXTRA.out.versions)
        ch_versions = ch_versions.mix(TABIX_DBSNP.out.versions)
        ch_versions = ch_versions.mix(BEDTOOLS_PAD_TARGET_BED.out.versions)
        ch_versions = ch_versions.mix(TABIX_BGZIPINDEX_PADDED_BED.out.versions)
        ch_versions = ch_versions.mix(GATK_BILT.out.versions)
        ch_versions = ch_versions.mix(GATK_ILT.out.versions)
        ch_versions = ch_versions.mix(CAT_CAT_BAIT.out.versions)
        ch_versions = ch_versions.mix(UNTAR_VEP_CACHE.out.versions)
        ch_versions = ch_versions.mix(RTGTOOLS_FORMAT.out.versions)

    emit:
        genome_bwa_index      = channel.empty().mix(ch_bwa, ch_sentieonbwa).collect()                        // channel: [ val(meta), path(index) ]
        genome_bwamem2_index  = ch_genome_bwamem2_index                                                      // channel: [ val(meta), path(index) ]
        genome_bwameme_index  = ch_genome_bwameme_index                                                      // channel: [ val(meta), path(index) ]
        genome_chrom_sizes    = GET_CHROM_SIZES.out.sizes.collect()                                          // channel: [ path(sizes) ]
        genome_fai            = ch_fai                                                                       // channel: [ val(meta), path(fai) ]
        genome_dict           = ch_genome_dict                                                               // channel: [ val(meta), path(dict) ]
        sdf                   = RTGTOOLS_FORMAT.out.sdf                                                      // channel: [ val (meta), path(intervals) ]
        mt_intervals          = ch_shiftfasta_mtintervals.intervals.collect()                                // channel: [ path(intervals) ]
        mt_bwa_index          = ch_mt_bwa_index                                                              // channel: [ val(meta), path(index) ]
        mt_bwamem2_index      = ch_mt_bwamem2_index                                                          // channel: [ val(meta), path(index) ]
        mt_dict               = ch_mt_dict                                                                   // channel: [ val(meta), path(dict) ]
        mt_fai                = ch_mt_fai                                                                    // channel: [ val(meta), path(fai) ]
        mt_fasta              = ch_mt_fasta                                                                  // channel: [ val(meta), path(fasta) ]
        mtshift_intervals     = ch_shiftfasta_mtintervals.shift_intervals.collect()                          // channel: [ path(intervals) ]
        mtshift_backchain     = GATK_SHIFTFASTA.out.shift_back_chain.collect()                               // channel: [ val(meta), path(backchain) ]
        mtshift_dict          = GATK_SHIFTFASTA.out.dict                                                     // channel: [ val(meta), path(dict) ]
        mtshift_fai           = GATK_SHIFTFASTA.out.shift_fai.collect()                                      // channel: [ val(meta), path(fai) ]
        mtshift_fasta         = GATK_SHIFTFASTA.out.shift_fa.collect()                                       // channel: [ val(meta), path(fasta) ]
        mtshift_bwa_index     = ch_mtshift_bwa_index                                                         // channel: [ val(meta), path(index) ]
        mtshift_bwamem2_index = ch_mtshift_bwamem2_index                                                     // channel: [ val(meta), path(index) ]
        gnomad_af_idx         = TABIX_GNOMAD_AF.out.tbi.collect()                                            // channel: [ val(meta), path(fasta) ]
        known_dbsnp_tbi       = TABIX_DBSNP.out.tbi.collect()                                                // channel: [ val(meta), path(fasta) ]
        target_bed            = ch_target_bed_gz_tbi.collect()                                               // channel: [ val(meta), path(bed), path(tbi) ]
        vcfanno_extra         = ch_vcfanno_extra.ifEmpty([[]])                                               // channel: [ [path(vcf), path(tbi)] ]
        bait_intervals        = CAT_CAT_BAIT.out.file_out.map{ _meta, inter -> inter}.collect().ifEmpty([[]])// channel: [ path(intervals) ]
        target_intervals      = GATK_BILT.out.interval_list.map{ _meta, inter -> inter}.collect()            // channel: [ path(interval_list) ]
        vep_resources         = UNTAR_VEP_CACHE.out.untar.map{ _meta, files -> [files]}.collect()            // channel: [ path(cache) ]
        versions              = ch_versions                                                                  // channel: [ path(versions.yml) ]

}
