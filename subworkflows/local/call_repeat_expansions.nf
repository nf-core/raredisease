//
// Run ExpansionHunter
//

params.expansionhunter_options = [:]

include { EXPANSIONHUNTER } from '../../modules/nf-core/modules/expansionhunter/main'  addParams( options: params.expansionhunter_options )

workflow CALL_REPEAT_EXPANSIONS {
    take:
        bam 				// channel: [ val(meta), path(bam), path(bai) ]
        fasta 			// channel: /path/to/genome.fasta
        variant_catalog		// channel: /path/to/variant_catalog.json

    main:

        ch_versions = Channel.empty()
        
        EXPANSIONHUNTER( bam, fasta, file(variant_catalog) )
        ch_versions = ch_versions.mix(EXPANSIONHUNTER.out.versions)

    emit:
        vcf         = EXPANSIONHUNTER.out.vcf       // channel: [ val(meta), path(*.vcf) ]
        versions    = ch_versions.ifEmpty(null)     // channel: [ versions.yml ]
}
