//
// Map to reference, fetch stats for each demultiplexed read pair, merge, mark duplicates, and index.
//

params.sambamba_options = [:]

include { SAMBAMBA } from '../../modules/nf-core/modules/sambamba/main'  addParams( options: params.sambamba_options )



workflow MITO {
    take:
        bam   // channel: [ val(meta), path(bam) ]
        

    main:
        // Map, sort, and index
        SAMBAMBA (bam )

        // Join the mapped bam + bai paths by their keys for stats
        // Get stats for each demultiplexed read pair.
        

    emit:
        sambamba_bam             = SAMBAMBA.out.bam         // channel: [ val(meta), [ marked_bam ] ]
        sambamba_bai             = SAMBAMBA.out.bai         // channel: [ val(meta), [ marked_bai ] ]

        // Collect versions
        sambamba_version       = SAMBAMBA.out.version      //      path: *.version.txt
}
