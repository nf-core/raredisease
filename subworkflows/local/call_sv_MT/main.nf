//
// Call SV MT
//

include { MT_DELETION  } from '../../../modules/local/mt_deletion_script'
include { MITOSALT     } from '../../../modules/local/mitosalt/main'
include { SEQTK_SAMPLE } from '../../../modules/nf-core/seqtk/sample/main'

workflow CALL_SV_MT {
    take:
        ch_reads        // channel: [mandatory] [ val(meta), [path(reads)] ]
        ch_bam_bai      // channel: [mandatory] [ val(meta), path(bam) ]
        ch_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_msconfig     // channel: [mandatory] [ path(msconfig) ]
        ch_genome       // channel: [mandatory] [ path(genome) ]

    main:
        ch_versions            = Channel.empty()
        ch_mitosalt_breakpoint = Channel.empty()
        ch_mitosalt_cluster    = Channel.empty()
        ch_subdepth            = params.mitosalt_depth
        ch_reads_subdepth      = ch_reads.map { meta, reads -> [meta, reads, ch_subdepth] }

        if (!(params.skip_tools && params.skip_tools.split(',').contains('mitosalt'))) {
            SEQTK_SAMPLE(ch_reads_subdepth)
            ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions.first())

            MITOSALT(SEQTK_SAMPLE.out.reads, ch_msconfig, ch_genome)
            ch_mitosalt_breakpoint = MITOSALT.out.breakpoint
            ch_mitosalt_cluster    = MITOSALT.out.cluster
            ch_versions            = ch_versions.mix(MITOSALT.out.versions.first())
        }
        MT_DELETION(ch_bam_bai, ch_fasta)

        ch_versions = ch_versions.mix(MT_DELETION.out.versions.first())

    emit:
        mitosalt_breakpoint = ch_mitosalt_breakpoint        // channel: [ val(meta), path(breakpoint) ]
        mitosalt_cluster    = ch_mitosalt_cluster           // channel: [ val(meta), path(cluster) ]
        mt_del_result       = MT_DELETION.out.mt_del_result // channel: [ val(meta), path(txt) ]
        versions            = ch_versions                   // channel: [ path(versions.yml) ]
}
