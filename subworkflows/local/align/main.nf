//
// Map to reference
//

include { ALIGN_BWA_BWAMEM2_BWAMEME                  } from '../align_bwa_bwamem2_bwameme'
include { ALIGN_MT                                   } from '../align_MT'
include { ALIGN_MT as ALIGN_MT_SHIFT                 } from '../align_MT'
include { ALIGN_SENTIEON                             } from '../align_sentieon'
include { CONVERT_MT_BAM_TO_FASTQ                    } from '../convert_mt_bam_to_fastq'
include { FASTP                                      } from '../../../modules/nf-core/fastp/main'
include { SAMTOOLS_VIEW as CONVERTTOBAM_CRAM         } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as CONVERTTOCRAM_ALTFILTERED } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as CONVERTTOCRAM_UNFILTERED  } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_EXCLUDE_ALT } from '../../../modules/nf-core/samtools/view/main'

workflow ALIGN {
    take:
        ch_alignments                 // channel: [optional] [ val(meta), [path(bam),path(bai)]  ]
        ch_genome_bwafastalignindex   // channel: [mandatory] [ val(meta), path(index) ]
        ch_genome_bwaindex            // channel: [mandatory] [ val(meta), path(index) ]
        ch_genome_bwamem2index        // channel: [mandatory] [ val(meta), path(index) ]
        ch_genome_bwamemeindex        // channel: [mandatory] [ val(meta), path(index) ]
        ch_genome_dictionary          // channel: [mandatory] [ val(meta), path(dict) ]
        ch_genome_fai                 // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_fasta               // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_input_reads                // channel: [optional] [ val(meta), [path(reads)]  ]
        ch_mt_bwaindex                // channel: [mandatory] [ val(meta), path(index) ]
        ch_mt_bwamem2index            // channel: [mandatory] [ val(meta), path(index) ]
        ch_mt_dictionary              // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mt_fai                     // channel: [mandatory] [ val(meta), path(fai) ]
        ch_mt_fasta                   // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_mtshift_bwaindex           // channel: [mandatory] [ val(meta), path(index) ]
        ch_mtshift_bwamem2index       // channel: [mandatory] [ val(meta), path(index) ]
        ch_mtshift_dictionary         // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mtshift_fai                // channel: [mandatory] [ val(meta), path(fai) ]
        ch_mtshift_fasta              // channel: [mandatory] [ val(meta), path(fasta) ]
        skip_fastp                    // boolean
        val_aligner                   //  string:  'bwa', 'bwafastalign', 'bwamem2', 'bwameme', or 'sentieon'
        val_analysis_type             //  string:  'wgs', 'wes', or 'mito'
        val_exclude_alt               // boolean
        val_extract_alignments        // boolean
        val_mt_aligner                //  string:  'bwa', 'bwamem2', or 'sentieon'
        val_platform                  //  string:  [mandatory] illumina or a different technology
        val_run_mt_for_wes            // boolean
        val_save_all_mapped_as_cram   // boolean
        val_save_noalt_mapped_as_cram // boolean

    main:
        ch_bwamem2_bam               = channel.empty()
        ch_bwamem2_bai               = channel.empty()
        ch_fastp_json                = channel.empty()
        ch_fastp_out                 = channel.empty()
        ch_genome_marked_cram        = channel.empty()
        ch_genome_marked_crai        = channel.empty()
        ch_markdup_metrics           = channel.empty()
        ch_mt_bam_bai                = channel.empty()
        ch_mt_bam_bai_gatksubwf      = channel.empty()
        ch_mtshift_bam_bai_gatksubwf = channel.empty()
        ch_sentieon_bam              = channel.empty()
        ch_sentieon_bai              = channel.empty()

        if (!skip_fastp) {
            FASTP (ch_input_reads.map {meta, reads -> return [meta, reads, []] }, false, false, false)
            ch_input_reads = FASTP.out.reads
            ch_fastp_json  = FASTP.out.json
            ch_fastp_out   = FASTP.out.json
                .mix(FASTP.out.html)
                .mix(FASTP.out.log)
                .mix(FASTP.out.reads.transpose())
                .mix(FASTP.out.reads_fail)
                .mix(FASTP.out.reads_merged)
        }

        //
        // If input is bam or cram
        //
        ch_alignments
            .map { meta, files ->
                def new_id   = meta.sample
                def new_meta = meta + [id:new_id, read_group:"\'@RG\\tID:" + new_id + "\\tPL:" + val_platform + "\\tSM:" + new_id + "\'"] - meta.subMap('lane')
                return [new_meta, files].flatten()
            }
            .branch { meta, _file1, _file2 ->
                bam:  meta.data_type == "bam"
                cram: meta.data_type == "cram"
            }
            .set { ch_alignments_by_type }

        ch_alignments_by_type.bam
            .multiMap { meta, bam, bai ->
                bam: [meta - meta.subMap('data_type'), bam]
                bai: [meta - meta.subMap('data_type'), bai]
            }
            .set { ch_input_aligned }

        ch_input_bam = ch_input_aligned.bam
        ch_input_bai = ch_input_aligned.bai

        CONVERTTOBAM_CRAM(
            ch_alignments_by_type.cram.map { meta, cram, crai -> [meta, cram, crai] },
            ch_genome_fasta.map { meta, fasta -> [meta, fasta, []] },
            [],
            'bai'
        )

        if (val_aligner.matches("bwamem2|bwa|bwameme|bwafastalign")) {
            ALIGN_BWA_BWAMEM2_BWAMEME (
                ch_genome_bwaindex,
                ch_genome_bwafastalignindex,
                ch_genome_bwamem2index,
                ch_genome_bwamemeindex,
                ch_genome_fai,
                ch_genome_fasta,
                ch_input_reads,
                val_aligner,
                val_extract_alignments,
                val_platform
            )
            ch_bwamem2_bam     = ALIGN_BWA_BWAMEM2_BWAMEME.out.marked_bam
            ch_bwamem2_bai     = ALIGN_BWA_BWAMEM2_BWAMEME.out.marked_bai
            ch_markdup_metrics = ALIGN_BWA_BWAMEM2_BWAMEME.out.markdup_metrics
        } else if (val_aligner.equals("sentieon")) {
            ALIGN_SENTIEON (
                ch_genome_bwaindex,
                ch_genome_fai,
                ch_genome_fasta,
                ch_input_reads,
                val_extract_alignments,
                val_platform
            )
            ch_sentieon_bam    = ALIGN_SENTIEON.out.marked_bam
            ch_sentieon_bai    = ALIGN_SENTIEON.out.marked_bai
            ch_markdup_metrics = ALIGN_SENTIEON.out.dedup_metrics
                .mix(ALIGN_SENTIEON.out.dedup_metrics_multiqc)
                .mix(ALIGN_SENTIEON.out.score)
        }

        ch_cram_converted_bam = CONVERTTOBAM_CRAM.out.bam.map { meta, bam -> [meta - meta.subMap('data_type'), bam] }
        ch_cram_converted_bai = CONVERTTOBAM_CRAM.out.bai.map { meta, bai -> [meta - meta.subMap('data_type'), bai] }

        ch_genome_marked_bam_initial     = channel.empty().mix(ch_bwamem2_bam, ch_sentieon_bam, ch_input_bam, ch_cram_converted_bam)
        ch_genome_marked_bai_initial     = channel.empty().mix(ch_bwamem2_bai, ch_sentieon_bai, ch_input_bai, ch_cram_converted_bai)
        ch_genome_marked_bam_bai_initial = ch_genome_marked_bam_initial.join(ch_genome_marked_bai_initial, failOnMismatch:true, failOnDuplicate:true)

        ch_branched = ch_genome_marked_bam_bai_initial.branch { meta, bam, bai ->
            exclude_alt: val_exclude_alt
                return [meta, bam, bai]
            keep_all: true
                return [meta, bam, bai]
        }

        SAMTOOLS_VIEW_EXCLUDE_ALT(
            ch_branched.exclude_alt,
            ch_genome_fasta.map { meta, fasta -> [meta, fasta, []] },
            [],
            'bai'
        )

        ch_genome_marked_bam     = val_exclude_alt ? SAMTOOLS_VIEW_EXCLUDE_ALT.out.bam : ch_branched.keep_all.map { meta, bam, _bai -> [meta, bam] }
        ch_genome_marked_bai     = val_exclude_alt ? SAMTOOLS_VIEW_EXCLUDE_ALT.out.bai : ch_branched.keep_all.map { meta, _bam, bai -> [meta, bai] }
        ch_genome_marked_bam_bai = ch_genome_marked_bam.join(ch_genome_marked_bai, failOnMismatch:true, failOnDuplicate:true)

        // PREPARING READS FOR MT ALIGNMENT

        if (val_analysis_type.matches("wgs|mito") || val_run_mt_for_wes) {
            CONVERT_MT_BAM_TO_FASTQ (
                ch_genome_marked_bam_bai,
                ch_genome_dictionary,
                ch_genome_fai,
                ch_genome_fasta
            )

            ALIGN_MT (
                ch_mt_bwaindex,
                ch_mt_bwamem2index,
                ch_mt_dictionary,
                ch_mt_fai,
                ch_mt_fasta,
                CONVERT_MT_BAM_TO_FASTQ.out.fastq,
                CONVERT_MT_BAM_TO_FASTQ.out.ubam,
                val_mt_aligner
            )

            ALIGN_MT_SHIFT (
                ch_mtshift_bwaindex,
                ch_mtshift_bwamem2index,
                ch_mtshift_dictionary,
                ch_mtshift_fai,
                ch_mtshift_fasta,
                CONVERT_MT_BAM_TO_FASTQ.out.fastq,
                CONVERT_MT_BAM_TO_FASTQ.out.ubam,
                val_mt_aligner
            )

            ch_mt_bam_bai                = CONVERT_MT_BAM_TO_FASTQ.out.bam_bai // Used for subsampling and SV calling
            ch_mt_bam_bai_gatksubwf      = ALIGN_MT.out.marked_bam
                                            .join(ALIGN_MT.out.marked_bai, failOnMismatch:true, failOnDuplicate:true) // Only for SNV calling
            ch_mtshift_bam_bai_gatksubwf = ALIGN_MT_SHIFT.out.marked_bam
                                            .join(ALIGN_MT_SHIFT.out.marked_bai, failOnMismatch:true, failOnDuplicate:true) // Only for SNV calling
        }

        if (val_save_noalt_mapped_as_cram) {
            CONVERTTOCRAM_ALTFILTERED( ch_genome_marked_bam_bai, ch_genome_fasta.map{meta, fasta -> return [meta, fasta, []]}, [], 'crai' )
            ch_genome_marked_cram = ch_genome_marked_cram.mix(CONVERTTOCRAM_ALTFILTERED.out.cram)
            ch_genome_marked_crai = ch_genome_marked_crai.mix(CONVERTTOCRAM_ALTFILTERED.out.crai)
        }

        if (val_save_all_mapped_as_cram) {
            CONVERTTOCRAM_UNFILTERED( ch_genome_marked_bam_bai_initial, ch_genome_fasta.map{meta, fasta -> return [meta, fasta, []]}, [], 'crai' )
            ch_genome_marked_cram = ch_genome_marked_cram.mix(CONVERTTOCRAM_UNFILTERED.out.cram)
            ch_genome_marked_crai = ch_genome_marked_crai.mix(CONVERTTOCRAM_UNFILTERED.out.crai)
        }

    emit:
        fastp_json                = ch_fastp_json                // channel: [ val(meta), path(json) ]
        fastp_out                 = ch_fastp_out                 // channel: [ val(meta), path(json|html|log|reads|reads_fail|reads_merged) ]
        genome_marked_cram        = ch_genome_marked_cram        // channel: [ val(meta), path(cram) ]
        genome_marked_crai        = ch_genome_marked_crai        // channel: [ val(meta), path(crai) ]
        genome_marked_bam         = ch_genome_marked_bam         // channel: [ val(meta), path(bam) ]
        genome_marked_bai         = ch_genome_marked_bai         // channel: [ val(meta), path(bai) ]
        genome_marked_bam_bai     = ch_genome_marked_bam_bai     // channel: [ val(meta), path(bam), path(bai) ]
        markdup_metrics           = ch_markdup_metrics           // channel: [ val(meta), path(metrics) ]
        mt_bam_bai                = ch_mt_bam_bai                // channel: [ val(meta), path(bam), path(bai) ]
        mt_bam_bai_gatksubwf      = ch_mt_bam_bai_gatksubwf      // channel: [ val(meta), path(bam), path(bai) ]
        mtshift_bam_bai_gatksubwf = ch_mtshift_bam_bai_gatksubwf // channel: [ val(meta), path(bam), path(bai) ]
}
