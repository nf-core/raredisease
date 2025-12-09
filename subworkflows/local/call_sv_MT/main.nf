//
// Call SV MT
//

include { MT_DELETION   } from '../../../modules/local/mt_deletion_script'
include { PREP_MITOSALT } from '../../../modules/local/prep_mitosalt/main'
include { MITOSALT      } from '../../../modules/local/mitosalt/main'
include { SEQTK_SAMPLE  } from '../../../modules/nf-core/seqtk/sample/main'

workflow CALL_SV_MT {
    take:
        ch_reads                     // channel: [mandatory] [ val(meta), [path(reads)] ]
        ch_bam_bai                   // channel: [mandatory] [ val(meta), path(bam) ]
        ch_genome_fasta              // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_hisat2index        // channel: [mandatory] [ val(meta), path(hisat2index) ]
	ch_genome_fai                // channel: [mandatory] [ val(meta), path(genomefai) ]
	ch_mt_lastdb                 // channel: [mandatory] [ val(meta), path(lastindex) ]
	ch_mt_fai                    // channel: [mandatory] [ val(meta), path(mtfai) ]
	ch_genome_chrsizes           // channel: [mandatory] [ val(meta), path(chrsizes) ]
	ch_mt_fasta                  // channel: [mandatory] [ val(meta), path(mtfasta) ]
        val_score_threshold          // string: [mandatory] mitosalt_score_threshold
        val_evalue_threshold         // string: [mandatory] mitosalt_evalue_threshold
        val_split_length             // string: [mandatory] mitosalt_split_length
        val_paired_distance          // string: [mandatory] mitosalt_paired_distance
        val_deletion_threshold_min   // string: [mandatory] mitosalt_del_threshold_min
        val_deletion_threshold_max   // string: [mandatory] mitosalt_del_threshold_max
        val_breakthreshold           // string: [mandatory] mitosalt_break_threshold
        val_cluster_threshold        // string: [mandatory] mitosalt_cluster_threshold
        val_breakspan                // string: [mandatory] mitosalt_break_span
        val_sizelimit                // string: [mandatory] mitosalt_size_limit
        val_hplimit                  // string: [mandatory] mitosalt_hp_limit
        val_flank                    // string: [mandatory] mitosalt_flank
        val_split_distance_threshold // string: [mandatory] mitosalt_split_dist_threshold
        ch_subdepth                  // channel: [mandatory] [ val(mitosalt_depth) ]
	ch_mito_name                 // channel: [mandatory] [ val(mito_name) ]
	val_exclude                  // string: [mandatory] mitosalt_exclude

    main:
        ch_versions = Channel.empty()

        if (!(params.skip_tools && params.skip_tools.split(',').contains('mitosalt'))) {
            ch_reads_subdepth      = ch_reads.concat(ch_subdepth).collect()
            SEQTK_SAMPLE (ch_reads_subdepth)
            ch_versions            = ch_versions.mix(SEQTK_SAMPLE.out.versions.first())

	    PREP_MITOSALT(
	        ch_genome_hisat2index,      
	        ch_genome_fai,              
	        ch_mt_lastdb,               
	        ch_mt_fai,                  
	        ch_genome_chrsizes,         
	        ch_mt_fasta,
		ch_mito_name,
		val_exclude,                
	        val_score_threshold,         
	        val_evalue_threshold,        
	        val_split_length,            
	        val_paired_distance,         
	        val_deletion_threshold_min,  
	        val_deletion_threshold_max,  
	        val_breakthreshold,          
	        val_cluster_threshold,       
	        val_breakspan,               
	        val_sizelimit,               
	        val_hplimit,                 
	        val_flank,                   
	        val_split_distance_threshold
	    ) 

            MITOSALT( 
		SEQTK_SAMPLE.out.reads, 
		PREP_MITOSALT.out.msconfig,
                ch_genome_hisat2index,
                ch_genome_fai,
                ch_mt_lastdb,
                ch_mt_fai,
                ch_genome_chrsizes,
                ch_mt_fasta
	    )
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
