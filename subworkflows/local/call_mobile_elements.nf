//
// A subworkflow to call mobile elements in the genome
//

include { RETROSEQ_CALL as RETROSEQ_CALL             } from '../../modules/local/retroseq/call/main'
include { RETROSEQ_DISCOVER as RETROSEQ_DISCOVER     } from '../../modules/local/retroseq/discover/main'
include { SAMTOOLS_VIEW as ME_SPLIT_ALIGNMENT        } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX as ME_INDEX_SPLIT_ALIGNMENT } from '../../modules/nf-core/samtools/index/main'

workflow CALL_MOBILE_ELEMENTS {

    take:
        ch_genome_bam_bai   // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_genome_fasta     // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_fai       // channel: [mandatory] [ val(meta), path(fai) ]
        ch_me_references    // channel: [mandatory] [path(tsv)]
        val_genome_build    // string: [mandatory] GRCh37 or GRCh38

    main:
        ch_versions = Channel.empty()

        // Building chromosome channel depending on genome version
        // TODO: Check how retroseq behaves when running chrY on female samples
        Channel.of(1..22, 'X', 'Y')
            .branch { it ->
                grch37: val_genome_build.equals('GRCh37')
                    return [it.toString()]
                grch38: val_genome_build.equals('GRCh38')
                    return ['chr' + it.toString()]
            }.set{ ch_chr_genome }
        ch_chr = ch_chr_genome.grch37.mix(ch_chr_genome.grch38)

        // Building one bam channel per chromosome and adding interval
        ch_genome_bam_bai
            .combine(ch_chr)
            .map {
                meta, bam, bai, chr ->
                return [ meta + [interval:chr], bam, bai ]
            }
            .set { ch_genome_bam_bai_interval }

        ME_SPLIT_ALIGNMENT(
            ch_genome_bam_bai_interval,
            [[:], []],
            []
        )

        ME_INDEX_SPLIT_ALIGNMENT( ME_SPLIT_ALIGNMENT.out.bam )

        ME_SPLIT_ALIGNMENT.out.bam
            .join(ME_INDEX_SPLIT_ALIGNMENT.out.bai, failOnMismatch:true)
            .set { ch_retroseq_input }

        ch_me_references
            .multiMap { type, path ->
                type: type
                path: path
            }
            .set { ch_me_reference_split }

        RETROSEQ_DISCOVER (
            ch_retroseq_input,
            ch_me_reference_split.path.collect(),
            ch_me_reference_split.type.collect()
        )

        RETROSEQ_DISCOVER.out.tab
            .join(ch_retroseq_input, failOnMismatch: true)
            .set { ch_retroseq_call_input }

        // Running retroseq call: clusters reads and checks on the breakpoints to decide whether a TEV is present
        RETROSEQ_CALL (
            ch_retroseq_call_input,
            ch_genome_fasta,
            ch_genome_fai
        )


        // Run vep to annotate

        // Run svdb to query against database

        // Filter and rank as done in findtroll? E.g. protein coding

        ch_versions = ch_versions.mix(ME_SPLIT_ALIGNMENT.out.versions).first()
        ch_versions = ch_versions.mix(ME_INDEX_SPLIT_ALIGNMENT.out.versions).first()
        ch_versions = ch_versions.mix(RETROSEQ_DISCOVER.out.versions).first()
        ch_versions = ch_versions.mix(RETROSEQ_CALL.out.versions).first()

    emit:
        versions = ch_versions              // channel: [ path(versions.yml) ]
}
