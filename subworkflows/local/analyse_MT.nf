//
// Analyse MT
//
include { CONVERT_MT_BAM_TO_FASTQ                      } from './mitochondria/convert_mt_bam_to_fastq'
include { ALIGN_AND_CALL_MT                            } from './mitochondria/align_and_call_MT'
include { ALIGN_AND_CALL_MT as ALIGN_AND_CALL_MT_SHIFT } from './mitochondria/align_and_call_MT'
include { PICARD_LIFTOVERVCF                           } from '../../modules/nf-core/picard/liftovervcf/main'
include { MERGE_ANNOTATE_MT                            } from './mitochondria/merge_annotate_MT'

workflow ANALYSE_MT {
    take:
        ch_bam                    // channel: [mandatory] [ val(meta), file(bam), file(bai) ]
        ch_cadd_header            // channel: [mandatory] [ path(txt) ]
        ch_cadd_scores            // channel: [mandatory] [ path(annotation) ]
        ch_genome_bwa_index       // channel: [mandatory] [ path(index) ]
        ch_genome_bwamem2_index   // channel: [mandatory] [ path(index) ]
        ch_genome_fasta_meta      // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fasta_no_meta   // channel: [mandatory] [ path(fasta) ]
        ch_genome_dict_meta       // channel: [mandatory] [ val(meta), path(dict) ]
        ch_genome_dict_no_meta    // channel: [mandatory] [ path(dict) ]
        ch_genome_fai             // channel: [mandatory] [ path(fai) ]
        ch_mt_intervals           // channel: [mandatory] [ path(interval_list) ]
        ch_shift_mt_bwa_index     // channel: [mandatory] [ path(index) ]
        ch_shift_mt_bwamem2_index // channel: [mandatory] [ path(index) ]
        ch_shift_mt_fasta         // channel: [mandatory] [ path(fasta) ]
        ch_shift_mt_dict          // channel: [mandatory] [ path(dict) ]
        ch_shift_mt_fai           // channel: [mandatory] [ path(fai) ]
        ch_shift_mt_intervals     // channel: [mandatory] [ path(interval_list) ]
        ch_shift_mt_backchain     // channel: [mandatory] [ path(back_chain) ]
        ch_vcfanno_resources      // channel: [mandatory] [ path(resources) ]
        ch_vcfanno_toml           // channel: [mandatory] [ path(toml) ]
        val_vep_genome            // string:  [mandatory] GRCh37 or GRCh38
        val_vep_cache_version     // string:  [mandatory] 107
        ch_vep_cache              // channel: [mandatory] [ path(cache) ]
        ch_case_info              // channel: [mandatory] [ val(case_info) ]

    main:
        ch_versions = Channel.empty()

        // PREPARING READS FOR MT ALIGNMENT
        CONVERT_MT_BAM_TO_FASTQ (
            ch_bam,
            ch_genome_fasta_meta,
            ch_genome_fai,
            ch_genome_dict_no_meta
        )

        // MT ALIGNMENT  AND VARIANT CALLING
        ALIGN_AND_CALL_MT (
            CONVERT_MT_BAM_TO_FASTQ.out.fastq,
            CONVERT_MT_BAM_TO_FASTQ.out.bam,
            ch_genome_bwa_index,
            ch_genome_bwamem2_index,
            ch_genome_fasta_no_meta,
            ch_genome_dict_no_meta,
            ch_genome_fai,
            ch_mt_intervals
        )

        ALIGN_AND_CALL_MT_SHIFT (
            CONVERT_MT_BAM_TO_FASTQ.out.fastq,
            CONVERT_MT_BAM_TO_FASTQ.out.bam,
            ch_shift_mt_bwa_index,
            ch_shift_mt_bwamem2_index,
            ch_shift_mt_fasta,
            ch_shift_mt_dict,
            ch_shift_mt_fai,
            ch_shift_mt_intervals
        )

        // LIFTOVER VCF FROM REFERENCE MT TO SHIFTED MT
        PICARD_LIFTOVERVCF (
            ALIGN_AND_CALL_MT_SHIFT.out.vcf,
            ch_genome_dict_no_meta,
            ch_shift_mt_backchain,
            ch_genome_fasta_no_meta
        )

        // MT MERGE AND ANNOTATE VARIANTS
        MERGE_ANNOTATE_MT(
            ALIGN_AND_CALL_MT.out.vcf,
            PICARD_LIFTOVERVCF.out.vcf_lifted,
            ch_cadd_header,
            ch_cadd_scores,
            ch_genome_fasta_no_meta,
            ch_genome_dict_meta,
            ch_genome_dict_no_meta,
            ch_genome_fai,
            ch_vcfanno_resources,
            ch_vcfanno_toml,
            val_vep_genome,
            val_vep_cache_version,
            ch_vep_cache,
            ch_case_info
        )

        ch_versions = ch_versions.mix(CONVERT_MT_BAM_TO_FASTQ.out.versions)
        ch_versions = ch_versions.mix(ALIGN_AND_CALL_MT.out.versions)
        ch_versions = ch_versions.mix(ALIGN_AND_CALL_MT_SHIFT.out.versions)
        ch_versions = ch_versions.mix(PICARD_LIFTOVERVCF.out.versions.first())
        ch_versions = ch_versions.mix(MERGE_ANNOTATE_MT.out.versions)

    emit:
        vcf           = MERGE_ANNOTATE_MT.out.vcf              // channel: [ val(meta), path(vcf) ]
        tbi           = MERGE_ANNOTATE_MT.out.tbi              // channel: [ val(meta), path(tbi) ]
        stats         = ALIGN_AND_CALL_MT.out.stats            // channel: [ val(meta), path(stats) ]
        filt_stats    = ALIGN_AND_CALL_MT.out.filt_stats       // channel: [ val(meta), path(tsv) ]
        stats_sh      = ALIGN_AND_CALL_MT_SHIFT.out.stats      // channel: [ val(meta), path(stats) ]
        filt_stats_sh = ALIGN_AND_CALL_MT_SHIFT.out.filt_stats // channel: [ val(meta), path(tsv) ]
        haplog        = MERGE_ANNOTATE_MT.out.haplog           // channel: [ val(meta), path(txt) ]
        report        = MERGE_ANNOTATE_MT.out.report           // channel: [ path(html) ]
        txt           = ALIGN_AND_CALL_MT.out.txt              // channel: [ val(meta), path(txt) ]
        html          = ALIGN_AND_CALL_MT.out.html             // channel: [ val(meta), path(html) ]
        txt_sh        = ALIGN_AND_CALL_MT_SHIFT.out.txt        // channel: [ val(meta), path(txt) ]
        html_sh       = ALIGN_AND_CALL_MT_SHIFT.out.html       // channel: [ val(meta), path(html) ]
        versions      = ch_versions                            // channel: [ path(versions.yml) ]
}
