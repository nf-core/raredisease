//
// Map to reference
//

include { FASTP                      } from '../../../modules/nf-core/fastp/main'
include { ALIGN_BWA_BWAMEM2_BWAMEME  } from '../align_bwa_bwamem2_bwameme'
include { ALIGN_SENTIEON             } from '../align_sentieon'
include { SAMTOOLS_VIEW              } from '../../../modules/nf-core/samtools/view/main'

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
        val_mbuffer_mem          // integer: [mandatory] memory in megabytes
        val_platform             // string:  [mandatory] illumina or a different technology
        val_sort_threads         // integer: [mandatory] number of sorting threads

    main:
        ch_bwamem2_bam        = Channel.empty()
        ch_bwamem2_bai        = Channel.empty()
        ch_fastp_json         = Channel.empty()
        ch_markdup_metrics    = Channel.empty()
        ch_versions           = Channel.empty()

        FASTP (ch_reads, [], false, false, false)
        ch_reads = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions)
        ch_fastp_json = FASTP.out.json

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
        }

        ch_genome_marked_bam = Channel.empty().mix(ch_bwamem2_bam, ch_input_bam)
        ch_genome_marked_bai = Channel.empty().mix(ch_bwamem2_bai, ch_input_bai)
        ch_genome_bam_bai    = ch_genome_marked_bam.join(ch_genome_marked_bai, failOnMismatch:true, failOnDuplicate:true)

    emit:
        fastp_json         = ch_fastp_json         // channel: [ val(meta), path(json) ]
        genome_marked_bam  = ch_genome_marked_bam  // channel: [ val(meta), path(bam) ]
        genome_marked_bai  = ch_genome_marked_bai  // channel: [ val(meta), path(bai) ]
        genome_bam_bai     = ch_genome_bam_bai     // channel: [ val(meta), path(bam), path(bai) ]
        markdup_metrics    = ch_markdup_metrics    // channel: [ val(meta), path(txt) ]
        versions           = ch_versions           // channel: [ path(versions.yml) ]
}
