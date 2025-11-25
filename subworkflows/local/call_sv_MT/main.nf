//
// Call SV MT
//

include { MT_DELETION  } from '../../../modules/local/mt_deletion_script'
include { MITOSALT     } from '../../../modules/local/mitosalt/main'
include { SEQTK_SAMPLE } from '../../../modules/nf-core/seqtk/sample/main'

workflow CALL_SV_MT {
    take:
        ch_reads              // channel: [mandatory] [ val(meta), [path(reads)] ]
        ch_bam_bai            // channel: [mandatory] [ val(meta), path(bam) ]
        ch_genome_fasta       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_hisat2index // channel: [mandatory] [ val(meta), path(hisat2index) ]
	ch_genome_fai         // channel: [mandatory] [ val(meta), path(genomefai) ]
	ch_mt_lastdb          // channel: [mandatory] [ val(meta), path(lastindex) ]
	ch_mt_fai             // channel: [mandatory] [ val(meta), path(mtfai) ]
	ch_genome_chrsizes    // channel: [mandatory] [ val(meta), path(chrsizes) ]
	ch_mt_fasta           // channel: [mandatory] [ val(meta), path(mtfasta) ]
        ch_msconfig           // channel: [mandatory] [ path(msconfig) ]
        ch_msref              // channel: [mandatory] [ path(genome) ]
        ch_subdepth           // channel: [mandatory] [ val(mitosalt_depth) ]

    main:
        ch_versions            = Channel.empty()

        if (!(params.skip_tools && params.skip_tools.split(',').contains('mitosalt'))) {
            ch_reads_subdepth      = ch_reads.map { meta, reads -> [meta, reads, '1000000'] }
            SEQTK_SAMPLE (ch_reads_subdepth)
            ch_versions            = ch_versions.mix(SEQTK_SAMPLE.out.versions.first())

//            MITOSALT(SEQTK_SAMPLE.out.reads, ch_msconfig, ch_msref, ch_genome_hisat2index, ch_genome_fai, ch_mt_lastdb, ch_mt_fai, ch_genome_chrsizes, ch_mt_fasta)
	  MITOSALT(SEQTK_SAMPLE.out.reads, ch_msconfig, ch_msref)  
          ch_versions            = ch_versions.mix(MITOSALT.out.versions.first())
        }
        MT_DELETION(ch_bam_bai, ch_genome_fasta)

        ch_versions = ch_versions.mix(MT_DELETION.out.versions.first())

    emit:
        mitosalt_breakpoint = MITOSALT.out.breakpoint       // channel: [ val(meta), path(breakpoint) ]
        mitosalt_cluster    = MITOSALT.out.cluster          // channel: [ val(meta), path(cluster) ]
        mt_del_result       = MT_DELETION.out.mt_del_result // channel: [ val(meta), path(txt) ]
        versions            = ch_versions                   // channel: [ path(versions.yml) ]
}
