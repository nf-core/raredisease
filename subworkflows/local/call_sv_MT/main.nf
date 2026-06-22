//
// Calls SV MT, concatenates FASTQs per sample, then runs MitoSalt and SaltShaker.
// Also detects the number of discordant pairs using the mitodel script.
//

include { CAT_FASTQ           } from '../../../modules/nf-core/cat/fastq/main'
include { MITOSALT            } from '../../../modules/local/mitosalt/main'
include { MT_DELETION         } from '../../../modules/local/mt_deletion_script'
include { PREP_MITOSALT       } from '../../../modules/local/prep_mitosalt/main'
include { SALTSHAKER_CALL     } from '../../../modules/nf-core/saltshaker/call/main'
include { SALTSHAKER_CLASSIFY } from '../../../modules/nf-core/saltshaker/classify/main'
include { FIND_CONCATENATE    } from '../../../modules/nf-core/find/concatenate/main'
include { SALTSHAKER_PLOT     } from '../../../modules/nf-core/saltshaker/plot/main'
include { SALTSHAKER_TO_HTML  } from '../../../modules/local/saltshaker_to_html/main'
include { SEQTK_SAMPLE        } from '../../../modules/nf-core/seqtk/sample/main'
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
        skip_mitosalt                         // Boolean
        val_heavy_strand_origin_start         // string: [mandatory] mitochondira_heavy_strand_origin_start
        val_heavy_strand_origin_end           // string: [mandatory] mitochondira_heavy_strand_origin_end
        val_light_strand_origin_start         // string: [mandatory] mitochondira_light_strand_origin_start
        val_light_strand_origin_end           // string: [mandatory] mitochondira_light_strand_origin_end
        val_mitochondria_length               // string: [mandatory] mito_length
        val_mitochondria_name                 // string: [mandatory] mito_name
        val_mitosalt_flank                    // string: [mandatory] mitosalt_flank
        val_mitosalt_heteroplasmy_limit       // string: [mandatory] mitosalt_heteroplasmy_limit

    main:
        ch_saltshaker_html  = channel.empty()
        ch_saltshaker_txt   = channel.empty()
        ch_saltshaker_vcf   = channel.empty()
        ch_saltshaker_plot  = channel.empty()

        if (!skip_mitosalt) {
            ch_reads
                .map { meta, reads ->
                        def sample_group_key = meta.sample
                        return [sample_group_key, meta, reads]
                }
                .groupTuple()
                .map { sample_id, meta_list, reads_list ->
                    def combined_meta = meta_list[0].clone()
                    combined_meta.id = sample_id
                    combined_meta.remove('lane')
                    combined_meta.remove('read_group')

                    def all_reads = reads_list.flatten()

                    [combined_meta, all_reads]
                }
                .set { ch_cat_fastq }

            CAT_FASTQ(ch_cat_fastq)

            ch_reads_subdepth = CAT_FASTQ.out.reads.combine(ch_subdepth)

            SEQTK_SAMPLE(ch_reads_subdepth)

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

            ch_prepmitosalt_config = PREP_MITOSALT.out.msconfig.collect()
            MITOSALT(
                SEQTK_SAMPLE.out.reads,
                ch_prepmitosalt_config,
                ch_genome_chrsizes,
                ch_genome_fai,
                ch_genome_hisat2index,
                ch_mt_fai,
                ch_mt_fasta,
                ch_mt_lastdb
            )

            MITOSALT.out.cluster
                .filter { _meta, out -> out.countLines() > 0 }
                .set { ch_cluster }

            MITOSALT.out.breakpoint
                .join(ch_cluster)
                .set { ch_saltshaker_in }

            // Saltshaker modules will only run if mitosalt called SVs and created a cluster file
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

            SALTSHAKER_CLASSIFY.out.txt
                .map { meta, txt ->
                    return [['id': meta.case_id], txt]
                    }
                .groupTuple()
                .set { ch_saltshaker_txts }

            FIND_CONCATENATE(
                ch_saltshaker_txts
            )
            ch_saltshaker_txt = FIND_CONCATENATE.out.file_out

            SALTSHAKER_TO_HTML(
                ch_saltshaker_txt
            )
            ch_saltshaker_html = SALTSHAKER_TO_HTML.out.classify_html

            SALTSHAKER_PLOT(
                SALTSHAKER_CLASSIFY.out.classify
            )
            ch_saltshaker_plot = SALTSHAKER_PLOT.out.plot

            // Only merge if Saltshaker produced VCFs; filter on the flat list before wrapping so
            // combine never sees an empty-spread issue (combining case_info with [[]] would produce
            // a 1-element channel item that breaks 2-param destructuring closures downstream).
            SALTSHAKER_CLASSIFY.out.vcf
                .map { _meta, vcf -> vcf }
                .filter { !it.isEmpty() } // .collect{} produces an ArrayBag which .filter can't be used on, so first convert to regular elements with filter, then collect those to filter the whole list
                .collect()
                .filter { !it.isEmpty() }
                .map { vcf_list -> [vcf_list] }
                .set { ch_vcf_file_list }

            ch_case_info
                .combine(ch_vcf_file_list)
                .set { ch_merge_input_vcfs }

            SVDB_MERGE( ch_merge_input_vcfs, [], true ).vcf
                .set {ch_saltshaker_vcf}
            // Saltshaker only runs if there are mitosalt calls. We update priority list when the
            // saltshaker vcf is created so the priority matches the list of vcfs that will be merged later
            ch_svcaller_priority = ch_svcaller_priority
                .concat(ch_saltshaker_vcf.map{ _it -> ["mitosalt"] })
                .collect()

        }
        MT_DELETION(ch_bam_bai, ch_genome_fasta)


    emit:
        saltshaker_vcf   = ch_saltshaker_vcf             // channel: [ val(meta), path(vcf) ]
        saltshaker_html  = ch_saltshaker_html            // channel: [ val(meta), path(html) ]
        saltshaker_plot  = ch_saltshaker_plot            // channel: [ val(meta), path(png) ]
        mt_del_result    = MT_DELETION.out.mt_del_result // channel: [ val(meta), path(txt) ]
        updated_priority = ch_svcaller_priority          // channel: [ val(["caller1", "caller2", ...]) ] - includes "mitosalt" if it ran
}
