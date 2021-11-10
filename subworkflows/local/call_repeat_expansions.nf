//
// Run ExpansionHunter
//

params.expansionhunter_options = [:]

include { EXPANSIONHUNTER } from '../../modules/nf-core/modules/expansionhunter/main'  addParams( options: params.expansionhunter_options )

workflow CALL_REPEAT_EXPANSIONS {
    take:
        bam 				// channel: [ val(meta), path(bam), path(bai) ]
        fasta 				// path: genome.fasta
        variant_catalog		// channel: /path/to/variant_catalog.json

    main:

        EXPANSIONHUNTER( bam, fasta, file(variant_catalog) )

    emit:
        vcf                       	= EXPANSIONHUNTER.out.vcf       // path: *.vcf
        expansionhunter_version     = EXPANSIONHUNTER.out.versions 	// path: versions.yml
}
