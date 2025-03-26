//
// Map to reference
//

include { FASTP                      } from '../../modules/nf-core/fastp/main'
include { ALIGN_BWA_BWAMEM2_BWAMEME  } from './align_bwa_bwamem2_bwameme'
include { ALIGN_SENTIEON             } from './align_sentieon'
include { SAMTOOLS_VIEW              } from '../../modules/nf-core/samtools/view/main'
include { ALIGN_MT                   } from './align_MT'
include { ALIGN_MT as ALIGN_MT_SHIFT } from './align_MT'
include { CONVERT_MT_BAM_TO_FASTQ    } from './convert_mt_bam_to_fastq'

workflow ALIGN {
    take:
        ch_reads                 // channel: [optional] [ val(meta), [path(reads)]  ]
        ch_alignments            // channel: [optional] [ val(meta), [path(bam),path(bai)]  ]
        ch_genome_fasta          // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai            // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_bwaindex       // channel: [mandatory] [ val(meta), path(index) ]
        ch_genome_bwamem2index   // channel: [mandatory] [ val(meta), path(index) ]
        ch_genome_bwamemeindex   // channel: [mandatory] [ val(meta), path(index) ]
        ch_genome_dictionary     // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mt_bwaindex           // channel: [mandatory] [ val(meta), path(index) ]
        ch_mt_bwamem2index       // channel: [mandatory] [ val(meta), path(index) ]
        ch_mt_dictionary         // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mt_fai                // channel: [mandatory] [ val(meta), path(fai) ]
        ch_mt_fasta              // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_mtshift_bwaindex      // channel: [mandatory] [ val(meta), path(index) ]
        ch_mtshift_bwamem2index  // channel: [mandatory] [ val(meta), path(index) ]
        ch_mtshift_dictionary    // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mtshift_fai           // channel: [mandatory] [ val(meta), path(fai) ]
        ch_mtshift_fasta         // channel: [mandatory] [ val(meta), path(fasta) ]
        val_mbuffer_mem          // integer: [mandatory] memory in megabytes
        val_platform             // string:  [mandatory] illumina or a different technology
        val_sort_threads         // integer: [mandatory] number of sorting threads

    main:
        ch_bwamem2_bam        = Channel.empty()
        ch_bwamem2_bai        = Channel.empty()
        ch_fastp_json         = Channel.empty()
        ch_markdup_metrics    = Channel.empty()
        ch_mt_bam_bai         = Channel.empty()
        ch_mt_marked_bam      = Channel.empty()
        ch_mt_marked_bai      = Channel.empty()
        ch_mtshift_bam_bai    = Channel.empty()
        ch_mtshift_marked_bam = Channel.empty()
        ch_mtshift_marked_bai = Channel.empty()
        ch_sentieon_bam       = Channel.empty()
        ch_sentieon_bai       = Channel.empty()
        ch_versions           = Channel.empty()

        if (!(params.skip_tools && params.skip_tools.split(',').contains('fastp'))) {
            FASTP (ch_reads, [], false, false, false)
            ch_reads = FASTP.out.reads
            ch_versions = ch_versions.mix(FASTP.out.versions)
            ch_fastp_json = FASTP.out.json
        }

        //
        // If input is bam
        //
        ch_alignments.map { meta, files ->
                    def new_id   = meta.sample
                    def new_meta = meta + [id:new_id, read_group:"\'@RG\\tID:" + new_id + "\\tPL:" + val_platform + "\\tSM:" + new_id + "\'"] - meta.subMap('lane')
                    return [new_meta, files].flatten()
                }
                .map { it -> [it[0], it[1]] }
                .set{ch_input_bam}

        ch_alignments.map { meta, files ->
                    def new_id   = meta.sample
                    def new_meta = meta + [id:new_id, read_group:"\'@RG\\tID:" + new_id + "\\tPL:" + val_platform + "\\tSM:" + new_id + "\'"] - meta.subMap('lane')
                    return [new_meta, files].flatten()
                }
                .map { it -> [it[0], it[2]] }
                .set{ch_input_bai}

        if (params.aligner.matches("bwamem2|bwa|bwameme")) {
            ALIGN_BWA_BWAMEM2_BWAMEME (             // Triggered when params.aligner is set as bwamem2 or bwa or bwameme
                ch_reads,
                ch_genome_bwaindex,
                ch_genome_bwamem2index,
                ch_genome_bwamemeindex,
                ch_genome_fasta,
                ch_genome_fai,
                val_mbuffer_mem,
                val_platform,
                val_sort_threads
            )
            ch_bwamem2_bam     = ALIGN_BWA_BWAMEM2_BWAMEME.out.marked_bam
            ch_bwamem2_bai     = ALIGN_BWA_BWAMEM2_BWAMEME.out.marked_bai
            ch_markdup_metrics = ALIGN_BWA_BWAMEM2_BWAMEME.out.metrics
            ch_versions        = ch_versions.mix(ALIGN_BWA_BWAMEM2_BWAMEME.out.versions)
        } else if (params.aligner.equals("sentieon")) {
            ALIGN_SENTIEON (                        // Triggered when params.aligner is set as sentieon
                ch_reads,
                ch_genome_fasta,
                ch_genome_fai,
                ch_genome_bwaindex,
                val_platform
            )
            ch_sentieon_bam    = ALIGN_SENTIEON.out.marked_bam
            ch_sentieon_bai    = ALIGN_SENTIEON.out.marked_bai
            ch_versions     = ch_versions.mix(ALIGN_SENTIEON.out.versions)
        }

        ch_genome_marked_bam = Channel.empty().mix(ch_bwamem2_bam, ch_sentieon_bam, ch_input_bam)
        ch_genome_marked_bai = Channel.empty().mix(ch_bwamem2_bai, ch_sentieon_bai, ch_input_bai)
        ch_genome_bam_bai    = ch_genome_marked_bam.join(ch_genome_marked_bai, failOnMismatch:true, failOnDuplicate:true)

        // PREPARING READS FOR MT ALIGNMENT

        if (params.analysis_type.matches("wgs|mito") || params.run_mt_for_wes) {
            CONVERT_MT_BAM_TO_FASTQ (
                ch_genome_bam_bai,
                ch_genome_fasta,
                ch_genome_fai,
                ch_genome_dictionary
            )

            ALIGN_MT (
                CONVERT_MT_BAM_TO_FASTQ.out.fastq,
                CONVERT_MT_BAM_TO_FASTQ.out.bam,
                ch_mt_bwaindex,
                ch_mt_bwamem2index,
                ch_mt_fasta,
                ch_mt_dictionary,
                ch_mt_fai
            )

            ALIGN_MT_SHIFT (
                CONVERT_MT_BAM_TO_FASTQ.out.fastq,
                CONVERT_MT_BAM_TO_FASTQ.out.bam,
                ch_mtshift_bwaindex,
                ch_mtshift_bwamem2index,
                ch_mtshift_fasta,
                ch_mtshift_dictionary,
                ch_mtshift_fai
            )

            ch_mt_marked_bam      = ALIGN_MT.out.marked_bam
            ch_mt_marked_bai      = ALIGN_MT.out.marked_bai
            ch_mt_bam_bai         = ch_mt_marked_bam.join(ch_mt_marked_bai, failOnMismatch:true, failOnDuplicate:true)
            ch_mtshift_marked_bam = ALIGN_MT_SHIFT.out.marked_bam
            ch_mtshift_marked_bai = ALIGN_MT_SHIFT.out.marked_bai
            ch_mtshift_bam_bai    = ch_mtshift_marked_bam.join(ch_mtshift_marked_bai, failOnMismatch:true, failOnDuplicate:true)
            ch_versions           = ch_versions.mix(ALIGN_MT.out.versions,
                                        ALIGN_MT_SHIFT.out.versions,
                                        CONVERT_MT_BAM_TO_FASTQ.out.versions)
        }

        if (params.save_mapped_as_cram) {
            SAMTOOLS_VIEW( ch_genome_bam_bai, ch_genome_fasta, [] )
            ch_versions   = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
        }
    emit:
        fastp_json         = ch_fastp_json         // channel: [ val(meta), path(json) ]
        genome_marked_bam  = ch_genome_marked_bam  // channel: [ val(meta), path(bam) ]
        genome_marked_bai  = ch_genome_marked_bai  // channel: [ val(meta), path(bai) ]
        genome_bam_bai     = ch_genome_bam_bai     // channel: [ val(meta), path(bam), path(bai) ]
        markdup_metrics    = ch_markdup_metrics    // channel: [ val(meta), path(txt) ]
        mt_marked_bam      = ch_mt_marked_bam      // channel: [ val(meta), path(bam) ]
        mt_marked_bai      = ch_mt_marked_bai      // channel: [ val(meta), path(bai) ]
        mt_bam_bai         = ch_mt_bam_bai         // channel: [ val(meta), path(bam), path(bai) ]
        mtshift_marked_bam = ch_mtshift_marked_bam // channel: [ val(meta), path(bam) ]
        mtshift_marked_bai = ch_mtshift_marked_bai // channel: [ val(meta), path(bai) ]
        mtshift_bam_bai    = ch_mtshift_bam_bai    // channel: [ val(meta), path(bam), path(bai) ]
        versions           = ch_versions           // channel: [ path(versions.yml) ]
}
