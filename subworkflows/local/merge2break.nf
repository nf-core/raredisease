//
// Merge bams and then break into contigs (e.g. chromosomes).
//

params.samtools_merge_options = [:]
params.samtools_view_options = [:]

include { SAMTOOLS_MERGE } from '../../modules/local/samtools/merge/main' addParams( options: params.samtools_merge_options )
include { SAMTOOLS_VIEW } from '../../modules/nf-core/modules/samtools/view/main' addParams( options: params.samtools_view_options )

workflow MERGE2BREAK {
    // chr = Channel.of(1..22, 'X', 'Y', 'M').map{ "chr" + it }
    chr = "chr20"
    // TODO: need to write a validation check for chromosomes present in bam files and then chop bam file for each chr present.

    take:
        fasta // channel: [mandatory] fasta
        bam // channel: [mandatory] [ val(meta), [ bam ] ]
        bai // channel: [mandatory for regional merge] [ val(meta), [ bai ] ]

    main:
        bam_bai = bam.join(bai)
        bam_bai.map{ meta, bam, bai ->
            new_meta = meta.clone()
            new_meta.id = new_meta.id.split('_')[0]
            [new_meta, bam, bai]
        }.groupTuple().branch{
            single: it[1].size() == 1
            multiple: it[1].size() > 1
        }.set{ bam_bwa }


        SAMTOOLS_MERGE ( bam_bwa.multiple, chr )
        // TODO: Singles? WELLL, we shall re-name instead of merge but still chop into contigs for each chr present.
        // SAMTOOLS_VIEW ( bam_bwa.single )

    emit:
        bam = SAMTOOLS_MERGE.out.merged_bam
        csi = SAMTOOLS_MERGE.out.csi
}
