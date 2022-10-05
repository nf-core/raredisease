//
// Run ExpansionHunter and Stranger
//

include { EXPANSIONHUNTER } from '../../modules/nf-core/expansionhunter/main'
include { STRANGER } from '../../modules/nf-core/stranger/main'

workflow CALL_REPEAT_EXPANSIONS {
    take:
        bam 			// channel: [ val(meta), path(bam), path(bai) ]
        fasta 			// channel: /path/to/genome.fasta
        variant_catalog	// channel: /path/to/variant_catalog.json

    main:
        ch_versions = Channel.empty()

        EXPANSIONHUNTER( bam, fasta, variant_catalog )
        ch_versions = ch_versions.mix(EXPANSIONHUNTER.out.versions)

        STRANGER ( EXPANSIONHUNTER.out.vcf, variant_catalog )
        ch_versions = ch_versions.mix(STRANGER.out.versions)

    emit:
        vcf         = STRANGER.out.vcf                    // channel: [ val(meta), path(*.vcf.gz) ]
        versions    = ch_versions.ifEmpty(null)     // channel: [ versions.yml ]
}
