//
// Map to reference
//

include { FASTP                      } from '../../modules/nf-core/fastp/main'
include { ALIGN_BWA_BWAMEM2          } from './alignment/align_bwa_bwamem2'
include { ALIGN_SENTIEON             } from './alignment/align_sentieon'
include { SAMTOOLS_VIEW              } from '../../modules/nf-core/samtools/view/main'
include { ALIGN_MT                   } from './alignment/align_MT'
include { ALIGN_MT as ALIGN_MT_SHIFT } from './alignment/align_MT'
include { CONVERT_MT_BAM_TO_FASTQ    } from './mitochondria/convert_mt_bam_to_fastq'

workflow ALIGN {
    take:
        ch_reads                 // channel: [mandatory] [ val(meta), [path(reads)]  ]
        ch_genome_fasta          // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai            // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_bwaindex       // channel: [mandatory] [ val(meta), path(index) ]
        ch_genome_bwamem2index   // channel: [mandatory] [ val(meta), path(index) ]
        ch_genome_dictionary     // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mtshift_bwaindex      // channel: [mandatory] [ val(meta), path(index) ]
        ch_mtshift_bwamem2index  // channel: [mandatory] [ val(meta), path(index) ]
        ch_mtshift_fasta         // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_mtshift_dictionary    // channel: [mandatory] [ val(meta), path(dict) ]
        ch_mtshift_fai           // channel: [mandatory] [ val(meta), path(fai) ]
        val_platform             // string:  [mandatory] illumina or a different technology

    main:
        ch_versions       = Channel.empty()
        ch_fastp_json     = Channel.empty()
        ch_bwamem2_bam    = Channel.empty()
        ch_sentieon_bam   = Channel.empty()
        ch_bwamem2_bai    = Channel.empty()
        ch_sentieon_bai   = Channel.empty()

        if (!params.skip_fastp) {
            FASTP (ch_reads, [], false, false)
            ch_reads = FASTP.out.reads
            ch_versions = ch_versions.mix(FASTP.out.versions)
            ch_fastp_json = FASTP.out.json
        }

        if (params.aligner.equals("bwamem2") || params.aligner.equals("bwa")) {
            ALIGN_BWA_BWAMEM2 (             // Triggered when params.aligner is set as bwamem2 or bwa
                ch_reads,
                ch_genome_bwaindex,
                ch_genome_bwamem2index,
                ch_genome_fasta,
                ch_genome_fai,
                val_platform
            )
            ch_bwamem2_bam = ALIGN_BWA_BWAMEM2.out.marked_bam
            ch_bwamem2_bai = ALIGN_BWA_BWAMEM2.out.marked_bai
            ch_versions   = ch_versions.mix(ALIGN_BWA_BWAMEM2.out.versions)
        } else if (params.aligner.equals("sentieon")) {
            ALIGN_SENTIEON (            // Triggered when params.aligner is set as sentieon
                ch_reads,
                ch_genome_fasta,
                ch_genome_fai,
                ch_genome_bwaindex,
                val_platform
            )
            ch_sentieon_bam = ALIGN_SENTIEON.out.marked_bam
            ch_sentieon_bai = ALIGN_SENTIEON.out.marked_bai
            ch_versions   = ch_versions.mix(ALIGN_SENTIEON.out.versions)
        }

        ch_genome_marked_bam = Channel.empty().mix(ch_bwamem2_bam, ch_sentieon_bam)
        ch_genome_marked_bai = Channel.empty().mix(ch_bwamem2_bai, ch_sentieon_bai)
        ch_genome_bam_bai    = ch_genome_marked_bam.join(ch_genome_marked_bai, failOnMismatch:true, failOnDuplicate:true)

        // PREPARING READS FOR MT ALIGNMENT
        CONVERT_MT_BAM_TO_FASTQ (
            ch_genome_bam_bai,
            ch_genome_fasta,
            ch_genome_fai,
            ch_genome_dictionary
        )

        ALIGN_MT (
            CONVERT_MT_BAM_TO_FASTQ.out.fastq,
            CONVERT_MT_BAM_TO_FASTQ.out.bam,
            ch_genome_bwaindex,
            ch_genome_bwamem2index,
            ch_genome_fasta,
            ch_genome_dictionary,
            ch_genome_fai
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

        ch_mt_marked_bam = ALIGN_MT.out.marked_bam
        ch_mt_marked_bai = ALIGN_MT.out.marked_bai
        ch_mt_bam_bai    = ch_mt_marked_bam.join(ch_mt_marked_bai, failOnMismatch:true, failOnDuplicate:true)

        ch_mtshift_marked_bam = ALIGN_MT_SHIFT.out.marked_bam
        ch_mtshift_marked_bai = ALIGN_MT_SHIFT.out.marked_bai
        ch_mtshift_bam_bai    = ch_mtshift_marked_bam.join(ch_mtshift_marked_bai, failOnMismatch:true, failOnDuplicate:true)

        if (params.save_mapped_as_cram) {
            SAMTOOLS_VIEW( ch_genome_bam_bai, ch_genome_fasta, [] )
            ch_versions   = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
        }
        ch_versions   = ch_versions.mix(ALIGN_MT.out.versions,
                                        ALIGN_MT_SHIFT.out.versions,
                                        CONVERT_MT_BAM_TO_FASTQ.out.versions)

    emit:
        fastp_json         = ch_fastp_json         // channel: [ val(meta), path(json) ]
        genome_marked_bam  = ch_genome_marked_bam  // channel: [ val(meta), path(bam) ]
        genome_marked_bai  = ch_genome_marked_bai  // channel: [ val(meta), path(bai) ]
        genome_bam_bai     = ch_genome_bam_bai     // channel: [ val(meta), path(bam), path(bai) ]
        mt_marked_bam      = ch_mt_marked_bam      // channel: [ val(meta), path(bam) ]
        mt_marked_bai      = ch_mt_marked_bai      // channel: [ val(meta), path(bai) ]
        mt_bam_bai         = ch_mt_bam_bai         // channel: [ val(meta), path(bam), path(bai) ]
        mtshift_marked_bam = ch_mtshift_marked_bam // channel: [ val(meta), path(bam) ]
        mtshift_marked_bai = ch_mtshift_marked_bai // channel: [ val(meta), path(bai) ]
        mtshift_bam_bai    = ch_mtshift_bam_bai    // channel: [ val(meta), path(bam), path(bai) ]
        versions           = ch_versions           // channel: [ path(versions.yml) ]
}
