//
// A subworkflow to annotate mobile elements in the genome
//

include { RETROSEQ_CALL as RETROSEQ_CALL } from '../../modules/local/retroseq/call/main'
include { RETROSEQ_DISCOVER as RETROSEQ_DISCOVER } from '../../modules/local/retroseq/discover/main'

workflow ANNOTATE_STRUCTURAL_VARIANTS {

    take:
        ch_bam              // channel: [mandatory] [ val(meta), path(bam) ]
        ch_bai              // channel: [mandatory] [ val(meta), path(bai) ]
        ch_bam_bai          // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_fasta     // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai       // channel: [mandatory] [ val(meta), path(fai) ]

    main:
        ch_versions = Channel.empty()

        // Running retroseq discover: identify discordant read pairs that might support a TE insertion
        RETROSEQ_DISCOVER (
            ch_bam_bai,
            "tab",
            "extrafiles"
        )

        // Running retroseq call: clusters reads and checks on the breakpoints to decide whether a TEV is present
        RETROSEQ_CALL (
            ch_bam_bai,
            ch_genome_fasta,
        )

        // Run vep to annotate

        // Run svdb to query against database

        // Filter and rank as done in findtroll? E.g. protein coding

        ch_versions = ch_versions.mix(RETROSEQ_CALL.out.versions)
        ch_versions = ch_versions.mix(RETROSEQ_DISCOVER.out.versions)

    emit:
        vcf_ann  = RETROSEQ_CALL.out.vcf // channel: [ val(meta), path(vcf) ]
        versions = ch_versions              // channel: [ path(versions.yml) ]
}

