
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
include { TABIX_BGZIPTABIX as TABIX_BGZIPINDEX_PADDED_BED    } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPINDEX_VCFANNOEXTRA  } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_VCFANNOEXTRA                  } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_DBSNP                         } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GNOMAD_AF                     } from '../../../modules/nf-core/tabix/tabix/main'
include { UNTAR as UNTAR_VEP_CACHE                           } from '../../../modules/nf-core/untar/main'

workflow PREPARE_REFERENCES {
    take:
        val_aligner                  // String: "bwa", "bwamem2", "sentieon" or "bwameme"
        val_analysis_type            // String: "wgs", "wes", or "mito"
        val_bwa                      // String: path to bwa index
        val_bwamem2                  // String: path to bwamem2 index
        val_bwameme                  // String: path to bwameme index
        val_fai                      // String: path to genome fasta index
        val_fasta                    // String: path to genome fasta
        val_gnomad_af                // String: path to gnomad allele frequency file
        val_gnomad_af_idx            // String: path to gnomad allele frequency file's index
        val_known_dbsnp              // String: path to dbsnp file
        val_known_dbsnp_tbi          // String: path to dbsnp file's index
        val_mtaligner                // String: "bwa", "bwamem2", or "sentieon"
        val_mtfasta                  // String: path to mitochondrial fasta
        val_run_mt_for_wes           // Boolean
        val_run_rtgvcfeval           // Boolean
        val_sdf                      // String: path to sdf file
        val_genome_dict              // String: path to genome dictionary
        val_target_bed               // String: path to target bed file
        val_vcfanno_extra            // String: path to additional annotation files used by vcfanno
        val_vep_cache                // String: path to vep cache folder


    main:
        ch_versions               = channel.empty()
        ch_bait_intervals         = channel.empty()
        ch_bwa                    = channel.empty()
        ch_genome_bwameme_index   = channel.empty()
        ch_genome_bwamem2_index   = channel.empty()
        ch_gnomad_af_idx          = channel.empty()
        ch_dbsnp                  = channel.value([[:],[]])
        ch_dbsnp_tbi              = channel.value([[:],[]])
        ch_mt_bwa_index           = channel.empty()
        ch_mt_bwamem2_index       = channel.empty()
        ch_mtshift_bwa_index      = channel.empty()
        ch_mtshift_bwamem2_index  = channel.empty()
        ch_mtshift_backchain      = channel.empty()
        ch_mtshift_dict           = channel.empty()
        ch_mtshift_fai            = channel.empty()
        ch_mtshift_fasta          = channel.empty()
        ch_mt_dict                = channel.empty()
        ch_mt_fai                 = channel.empty()
        ch_mt_fasta               = channel.empty()
        ch_sdf                    = channel.empty()
        ch_shiftfasta_mtintervals = channel.empty()
        ch_target_bed_gz_tbi      = channel.value([[:],[],[]])
        ch_target_intervals       = channel.empty()
        ch_vcfanno_extra          = channel.value([[]])
        ch_vep_resources          = channel.value([[]])

        ch_genome_fasta = channel.fromPath(val_fasta).map { it -> [[id:it.simpleName], it] }.collect()
        //
        // Genome indices
        //
        if (!val_fai) {
            ch_genome_fai = SAMTOOLS_FAIDX_GENOME(ch_genome_fasta, [[:],[]], false).fai.collect()
        } else {
            ch_genome_fai = channel.fromPath(val_fai).map {it -> [[id:it.simpleName], it]}.collect()
        }

        GET_CHROM_SIZES( ch_genome_fai )

        if (!val_genome_dict) {
            ch_genome_dict = GATK_SD(ch_genome_fasta).dict.collect()
            ch_versions    = ch_versions.mix(GATK_SD.out.versions)
        } else {
            ch_genome_dict = channel.fromPath(val_genome_dict).map {it -> [[id:it.simpleName], it]}.collect()
        }
        //
        // Genome alignment indices
        //
        if (!val_bwa) {
            if (!val_aligner.equals("sentieon") || val_mtaligner.equals("bwa")) {
                ch_bwa      = BWA_INDEX_GENOME(ch_genome_fasta).index.collect()
                ch_versions = ch_versions.mix(BWA_INDEX_GENOME.out.versions)
            }
            if (val_aligner.equals("sentieon") || val_mtaligner.equals("sentieon")) {
                ch_bwa      = SENTIEON_BWAINDEX_GENOME(ch_genome_fasta).index.collect()
            }
        } else if (val_bwa) {
            ch_bwa = channel.fromPath(val_bwa).map {it -> [[id:it.simpleName], it]}.collect()
        }

        if (!val_bwamem2 && (val_aligner.equals("bwamem2") || val_mtaligner.equals("bwamem2"))) {
            ch_genome_bwamem2_index = BWAMEM2_INDEX_GENOME(ch_genome_fasta).index.collect()
            ch_versions             = ch_versions.mix(BWAMEM2_INDEX_GENOME.out.versions)
        } else if (val_bwamem2) {
            ch_genome_bwamem2_index = channel.fromPath(val_bwamem2).map {it -> [[id:it.simpleName], it]}.collect()
        }

        if (!val_bwamem2 && val_aligner.equals("bwameme")) {
            ch_genome_bwameme_index = BWAMEME_INDEX_GENOME(ch_genome_fasta).index.collect()
            ch_versions             = ch_versions.mix(BWAMEME_INDEX_GENOME.out.versions)
        } else if (val_bwameme) {
            ch_genome_bwameme_index = channel.fromPath(val_bwameme).map {it -> [[id:it.simpleName], it]}.collect()
        }
        //
        // MT genome indices
        //
        if (val_analysis_type.matches("wgs|mito") || val_run_mt_for_wes) {
            if (!val_mtfasta) {
                ch_mt_fasta = SAMTOOLS_EXTRACT_MT(ch_genome_fasta, ch_genome_fai, false).fa.collect()
            } else {
                ch_mt_fasta = channel.fromPath(val_mtfasta).map { it -> [[id:it.simpleName], it] }.collect()
            }

            ch_mt_fai  = SAMTOOLS_FAIDX_MT(ch_mt_fasta, [[:],[]], false).fai.collect()
            ch_mt_dict = GATK_SD_MT(ch_mt_fasta).dict.collect()

            GATK_SHIFTFASTA(ch_mt_fasta, ch_mt_fai, ch_mt_dict)

            ch_mtshift_backchain     = GATK_SHIFTFASTA.out.shift_back_chain.collect()
            ch_mtshift_dict          = GATK_SHIFTFASTA.out.dict.collect()
            ch_mtshift_fai           = GATK_SHIFTFASTA.out.shift_fai.collect()
            ch_mtshift_fasta         = GATK_SHIFTFASTA.out.shift_fa.collect()

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

            ch_versions = ch_versions.mix(GATK_SD_MT.out.versions,
                                            GATK_SHIFTFASTA.out.versions)
        }
        //
        // MT alignment indices
        //
        if ((val_analysis_type.matches("wgs|mito") || val_run_mt_for_wes) && val_mtaligner.equals("bwamem2")) {
            ch_mt_bwamem2_index      = BWAMEM2_INDEX_MT(ch_mt_fasta).index.collect()
            ch_mtshift_bwamem2_index = BWAMEM2_INDEX_MT_SHIFT(ch_mtshift_fasta).index.collect()
            ch_versions              = ch_versions.mix(BWAMEM2_INDEX_MT.out.versions, BWAMEM2_INDEX_MT_SHIFT.out.versions)
        }
        if ((val_analysis_type.matches("wgs|mito") || val_run_mt_for_wes) && val_mtaligner.equals("bwa")) {
            ch_mt_bwa_index          = BWA_INDEX_MT(ch_mt_fasta).index.collect()
            ch_mtshift_bwa_index     = BWA_INDEX_MT_SHIFT(ch_mtshift_fasta).index.collect()
            ch_versions              = ch_versions.mix(BWA_INDEX_MT.out.versions, BWA_INDEX_MT_SHIFT.out.versions)
        }
        if ((val_analysis_type.matches("wgs|mito") || val_run_mt_for_wes) && val_mtaligner.equals("sentieon")) {
            ch_mt_bwa_index          = SENTIEON_BWAINDEX_MT(ch_mt_fasta).index.collect()
            ch_mtshift_bwa_index     = SENTIEON_BWAINDEX_MT_SHIFT(ch_mtshift_fasta).index.collect()
        }
        //
        // Vcf, tab and bed indices
        //
        if (val_known_dbsnp) {
            ch_dbsnp = channel.fromPath(val_known_dbsnp).map{ it -> [[id:it.simpleName], it] }.collect()
            if (val_known_dbsnp_tbi) {
                ch_dbsnp_tbi = channel.fromPath(val_known_dbsnp_tbi).map{ it -> [[id:it.simpleName], it] }.collect()
            } else {
                ch_dbsnp_tbi = TABIX_DBSNP(ch_dbsnp).index.collect()
            }
        }

        if (val_gnomad_af) {
            ch_gnomad_af = channel.fromPath(val_gnomad_af).map{ it -> [[id:it.simpleName], it] }.collect()
            if (val_gnomad_af_idx) {
                ch_gnomad_idx = channel.fromPath(val_gnomad_af_idx).map{ it -> [[id:it.simpleName], it] }.collect()
            } else {
                ch_gnomad_idx = TABIX_GNOMAD_AF(ch_gnomad_af).index.collect()
            }
            ch_gnomad_af_idx = ch_gnomad_af.join(ch_gnomad_idx).map {_meta, tab, idx -> [tab,idx]}.collect()
        }
        //
        // Index target bed file
        //
        if (val_target_bed) {
            ch_target_bed = channel.fromPath(val_target_bed).map{ it -> [[id:it.simpleName], it] }.collect()

            BEDTOOLS_PAD_TARGET_BED(
                ch_target_bed,
                ch_genome_fai.map { _meta, fai -> return fai }
            )
            ch_target_bed_gz_tbi = TABIX_BGZIPINDEX_PADDED_BED(BEDTOOLS_PAD_TARGET_BED.out.bed).gz_index

            ch_target_intervals = GATK_BILT(ch_target_bed, ch_genome_dict).interval_list.map{ _meta, inter -> inter}.collect()

            GATK_ILT(GATK_BILT.out.interval_list)

            GATK_ILT.out.interval_list
                .collect{ _meta, list -> list }
                .map { list ->
                    def meta = list.toString().split("_split")[0].split("/")[-1] + "_bait.intervals_list"
                    return [[id:meta], list]
                }
                .set { ch_bait_intervals_cat_in }

            ch_bait_intervals = CAT_CAT_BAIT ( ch_bait_intervals_cat_in ).file_out.map {_meta, inter -> inter}.collect().ifEmpty([[]])

            ch_versions = ch_versions.mix(GATK_BILT.out.versions, GATK_ILT.out.versions)
        }
        //
        // Prepare vcfanno extra files
        //
        if (val_vcfanno_extra) {
            ch_vcfanno_tabix_in = channel.fromPath(val_vcfanno_extra).map { it -> [[id:it.baseName], it] }

            if (val_vcfanno_extra.endsWith(".gz")) {
                TABIX_VCFANNOEXTRA(ch_vcfanno_tabix_in).index
                    .join(ch_vcfanno_tabix_in)
                    .map { _meta, tbi, vcf -> return [[vcf,tbi]]}
                    .set {ch_vcfanno_extra}
            } else {
                TABIX_BGZIPINDEX_VCFANNOEXTRA(ch_vcfanno_tabix_in)
                channel.empty()
                    .mix(TABIX_BGZIPINDEX_VCFANNOEXTRA.out.gz_index, TABIX_BGZIPINDEX_VCFANNOEXTRA.out.gz_csi)
                    .map { _meta, vcf, index -> return [[vcf,index]] }
                    .set {ch_vcfanno_extra}
            }
        }
        //
        // Prepare vep cache files
        //
        if (val_vep_cache) {
            if (val_vep_cache.endsWith("tar.gz")) {
                ch_vep_resources = UNTAR_VEP_CACHE (channel.fromPath(val_vep_cache).map { it -> [[id:'vep_cache'], it] }.collect()).untar.map{ _meta, files -> [files]}.collect()
                ch_versions      = ch_versions.mix(UNTAR_VEP_CACHE.out.versions)
            } else {
                ch_vep_resources = channel.fromPath(val_vep_cache).collect()
            }
        }
        //
        // RTG tools
        //
        if (!val_sdf && val_run_rtgvcfeval) {
            ch_genome_fasta.map { meta, fasta -> return [meta, fasta, [], [] ] }
                .set {ch_rtgformat_in}
            ch_sdf      = RTGTOOLS_FORMAT(ch_rtgformat_in).out.sdf
            ch_versions = ch_versions.mix(RTGTOOLS_FORMAT.out.versions)
        }

        // Gather versions
        ch_versions = ch_versions.mix(GET_CHROM_SIZES.out.versions)

    emit:
        bait_intervals        = ch_bait_intervals                                   // channel:[ path(intervals) ]
        dbsnp                 = ch_dbsnp                                            // channel:[ val(meta), path(dbsnp) ]
        dbsnp_tbi             = ch_dbsnp_tbi                                        // channel:[ val(meta), path(dbsnp_idx) ]
        genome_bwa_index      = ch_bwa                                              // channel:[ val(meta), path(index) ]
        genome_bwamem2_index  = ch_genome_bwamem2_index                             // channel:[ val(meta), path(index) ]
        genome_bwameme_index  = ch_genome_bwameme_index                             // channel:[ val(meta), path(index) ]
        genome_chrom_sizes    = GET_CHROM_SIZES.out.sizes.collect()                 // channel:[ path(sizes) ]
        genome_fai            = ch_genome_fai                                       // channel:[ val(meta), path(fai) ]
        genome_fasta          = ch_genome_fasta                                     // channel:[ val(meta), path(fasta) ]
        genome_dict           = ch_genome_dict                                      // channel:[ val(meta), path(dict) ]
        gnomad_af_idx         = ch_gnomad_af_idx                                    // channel:[ val(gnomad), path(idx) ]
        mt_bwa_index          = ch_mt_bwa_index                                     // channel:[ val(meta), path(index) ]
        mt_bwamem2_index      = ch_mt_bwamem2_index                                 // channel:[ val(meta), path(index) ]
        mt_dict               = ch_mt_dict                                          // channel:[ val(meta), path(dict) ]
        mt_fai                = ch_mt_fai                                           // channel:[ val(meta), path(fai) ]
        mt_fasta              = ch_mt_fasta                                         // channel:[ val(meta), path(fasta) ]
        mt_intervals          = ch_shiftfasta_mtintervals.intervals.collect()       // channel:[ path(intervals) ]
        mtshift_backchain     = ch_mtshift_backchain                                // channel:[ val(meta), path(backchain) ]
        mtshift_bwa_index     = ch_mtshift_bwa_index                                // channel:[ val(meta), path(index) ]
        mtshift_bwamem2_index = ch_mtshift_bwamem2_index                            // channel:[ val(meta), path(index) ]
        mtshift_dict          = ch_mtshift_dict                                     // channel:[ val(meta), path(dict) ]
        mtshift_fai           = ch_mtshift_fai                                      // channel:[ val(meta), path(fai) ]
        mtshift_fasta         = ch_mtshift_fasta                                    // channel:[ val(meta), path(fasta) ]
        mtshift_intervals     = ch_shiftfasta_mtintervals.shift_intervals.collect() // channel:[ path(intervals) ]
        sdf                   = ch_sdf                                              // channel:[ val (meta), path(sdf) ]
        target_bed            = ch_target_bed_gz_tbi.collect()                      // channel:[ val(meta), path(bed), path(tbi) ]
        target_intervals      = ch_target_intervals                                 // channel:[ path(interval_list) ]
        vcfanno_extra         = ch_vcfanno_extra                                    // channel:[ [path(vcf), path(tbi)] ]
        vep_resources         = ch_vep_resources                                    // channel:[ path(cache) ]
        versions              = ch_versions                                         // channel:[ path(versions.yml) ]

}
