//
// A quality check subworkflow for processed bams.
//

params.picard_collectmultiplemetrics_options = [:]
params.qualimap_bamqc_options = [:]

include { PICARD_COLLECTMULTIPLEMETRICS } from '../../modules/nf-core/modules/picard/collectmultiplemetrics/main'  addParams( options: params.picard_collectmultiplemetrics_options )
include { QUALIMAP_BAMQC } from '../../modules/nf-core/modules/qualimap/bamqc/main'  addParams( options: params.qualimap_bamqc_options )

workflow QC_BAM {

    take:
        bam     // channel: [ val(meta), path(bam) ]
        fasta   // path: genome.fasta
        gff     // path: file.gff
        use_gff // boolean

    main:
        PICARD_COLLECTMULTIPLEMETRICS ( bam, fasta )
        QUALIMAP_BAMQC ( bam, gff, use_gff )

    emit:
        multiple_metrics        = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
        qualimap_results        = QUALIMAP_BAMQC.out.results
}
