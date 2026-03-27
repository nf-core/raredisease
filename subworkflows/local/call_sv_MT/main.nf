//
// Call SV MT
//

include { MT_DELETION         } from '../../../modules/local/mt_deletion_script'
include { PREP_MITOSALT       } from '../../../modules/local/prep_mitosalt/main'
include { MITOSALT            } from '../../../modules/local/mitosalt/main'
include { SEQTK_SAMPLE        } from '../../../modules/nf-core/seqtk/sample/main'
include { SALTSHAKER_CALL     } from '../../../modules/nf-core/saltshaker/call/main'
include { SALTSHAKER_CLASSIFY } from '../../../modules/nf-core/saltshaker/classify/main'
include { SALTSHAKER_PLOT     } from '../../../modules/nf-core/saltshaker/plot/main'
include { SVDB_MERGE          } from '../../../modules/nf-core/svdb/merge/main'

workflow CALL_SV_MT {
    take:
        ch_bam_bai                            // channel: [mandatory] [ val(meta), path(bam) ]
        ch_case_info                          // channel: [mandatory] [ val(case_info) ]
        ch_genome_chrsizes                    // channel: [mandatory] [ path(chrsizes) ]
        ch_genome_fai                         // channel: [mandatory] [ val(meta), path(genomefai) ]
        ch_genome_fasta                       // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_hisat2index                 // channel: [mandatory] [ val(meta), path(hisat2index) ]
        ch_mt_fai                             // channel: [mandatory] [ val(meta), path(mtfai) ]
        ch_mt_fasta                           // channel: [mandatory] [ val(meta), path(mtfasta) ]
        ch_mt_lastdb                          // channel: [mandatory] [ val(meta), path(lastindex) ]
        ch_reads                              // channel: [mandatory] [ val(meta), [path(reads)] ]
        ch_subdepth                           // channel: [mandatory] [ val(mitosalt_depth) ]
        ch_svcaller_priority                  // channel: [mandatory] [ val(["var caller tag 1", ...]) ]
        ch_mitosalt_config                    // channel: [mandatory] [val(mitosalt_breakspan),val(mitosalt_breakthreshold),...,val(mitosalt_split_length)]
        val_heavy_strand_origin_start         // string: [mandatory] mitochondira_heavy_strand_origin_start
        val_heavy_strand_origin_end           // string: [mandatory] mitochondira_heavy_strand_origin_end
        val_light_strand_origin_start         // string: [mandatory] mitochondira_light_strand_origin_start
        val_light_strand_origin_end           // string: [mandatory] mitochondira_light_strand_origin_end
        val_mitochondria_length               // string: [mandatory] mito_length
        val_mitochondria_name                 // string: [mandatory] mito_name
        val_mitosalt_flank                    // string: [mandatory] mitosalt_flank
        val_mitosalt_heteroplasmy_limit       // string: [mandatory] mitosalt_heteroplasmy_limit

    main:
        ch_saltshaker_txt   = channel.empty()
        ch_saltshaker_vcf   = channel.empty()
        ch_saltshaker_plot  = channel.empty()

        if (!(params.skip_tools && params.skip_tools.split(',').contains('mitosalt'))) {
            ch_reads_subdepth      = ch_reads.combine(ch_subdepth)

            SEQTK_SAMPLE (ch_reads_subdepth)

            PREP_MITOSALT(
                ch_genome_chrsizes,
                ch_genome_fai,
                ch_genome_hisat2index,
                ch_mitosalt_config,
                ch_mt_fai,
                ch_mt_fasta,
                ch_mt_lastdb,
                val_mitosalt_flank,
                val_mitosalt_heteroplasmy_limit,
                val_mitochondria_name
            )

            MITOSALT(
                SEQTK_SAMPLE.out.reads,
                PREP_MITOSALT.out.msconfig,
                ch_genome_chrsizes,
                ch_genome_fai,
                ch_genome_hisat2index,
                ch_mt_fai,
                ch_mt_fasta,
                ch_mt_lastdb
            )

            MITOSALT.out.cluster
                .filter{ _meta, out -> out.countLines() > 0 }
                .set{ ch_cluster }

            MITOSALT.out.breakpoint
                .join(ch_cluster)
                .set{ ch_saltshaker_in }

            SALTSHAKER_CALL(
                ch_saltshaker_in,
                ch_mt_fasta,
                val_mitosalt_flank,
                val_mitosalt_heteroplasmy_limit,
                val_mitochondria_length,
                val_heavy_strand_origin_start,
                val_heavy_strand_origin_end,
                val_light_strand_origin_start,
                val_light_strand_origin_end
            )

            SALTSHAKER_CLASSIFY(
                SALTSHAKER_CALL.out.call,
                val_mitochondria_name
            )
            ch_saltshaker_txt = SALTSHAKER_CLASSIFY.out.txt
            ch_saltshaker_vcf = SALTSHAKER_CLASSIFY.out.vcf

            SALTSHAKER_PLOT(
                SALTSHAKER_CLASSIFY.out.classify
            )
            ch_saltshaker_plot = SALTSHAKER_PLOT.out.plot

            SALTSHAKER_CLASSIFY.out.vcf
                .collect{ _meta, vcf -> vcf }
                .toList()
                .set { ch_vcf_file_list }

            // Only set the channel if ch_vcf_file_list was made (ie combined channel has two elements)
            ch_case_info
                .combine(ch_vcf_file_list)
                .filter{ it -> it.size() == 2}
                .set { ch_merge_input_vcfs }

            SVDB_MERGE ( ch_merge_input_vcfs, [], true ).vcf
                .set {ch_saltshaker_vcf}
            // Update priority list when we know saltshaker will run (ie saltshaker vcf is created)
            // Updated priority list will be used when saltshaker vcf is merged with other SV vcfs
            ch_svcaller_priority = ch_svcaller_priority
                .concat(ch_saltshaker_vcf.map{ _meta -> ["mitosalt"] })
                .collect()

        }
        MT_DELETION(ch_bam_bai, ch_genome_fasta)

        ch_publish = ch_saltshaker_txt
                        .mix(ch_saltshaker_vcf)
                        .mix(ch_saltshaker_plot)
                        .mix(MT_DELETION.out.mt_del_result)
                        .map { meta, value -> ['call_sv/', [meta, value]] }

    emit:
        saltshaker_txt   = ch_saltshaker_txt             // channel: [ val(meta), path(txt) ]
        saltshaker_vcf   = ch_saltshaker_vcf             // channel: [ val(meta), path(vcf) ]
        saltshaker_plot  = ch_saltshaker_plot            // channel: [ val(meta), path(png) ]
        mt_del_result    = MT_DELETION.out.mt_del_result // channel: [ val(meta), path(txt) ]
        updated_priority = ch_svcaller_priority          // channel: [ val(["caller1", "caller2", ...]) ] - includes "mitosalt" if it ran
        publish          = ch_publish                    // channel: [ val(destination), val(value) ]
}
