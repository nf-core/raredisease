//
// Call SV MT
//

include { MT_DELETION   } from '../../../modules/local/mt_deletion_script'
include { PREP_MITOSALT } from '../../../modules/local/prep_mitosalt/main'
include { MITOSALT      } from '../../../modules/local/mitosalt/main'
include { SEQTK_SAMPLE  } from '../../../modules/nf-core/seqtk/sample/main'

workflow CALL_SV_MT {
    take:
        ch_reads                    // channel: [mandatory] [ val(meta), [path(reads)] ]
        ch_bam_bai                  // channel: [mandatory] [ val(meta), path(bam) ]
        ch_genome_fasta             // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_genome_hisat2index       // channel: [mandatory] [ val(meta), path(hisat2index) ]
	ch_genome_fai               // channel: [mandatory] [ val(meta), path(genomefai) ]
	ch_mt_lastdb                // channel: [mandatory] [ val(meta), path(lastindex) ]
	ch_mt_fai                   // channel: [mandatory] [ val(meta), path(mtfai) ]
	ch_genome_chrsizes          // channel: [mandatory] [ val(meta), path(chrsizes) ]
	ch_mt_fasta                 // channel: [mandatory] [ val(meta), path(mtfasta) ]
        ch_score_threshold          // channel: [mandatory] [ val(mitosalt_score_threshold) ]
        ch_evalue_threshold         // channel: [mandatory] [ val(mitosalt_evalue_threshold) ]
        ch_split_length             // channel: [mandatory] [ val(mitosalt_split_length) ]
        ch_paired_distance          // channel: [mandatory] [ val(mitosalt_paired_distance) ]
        ch_deletion_threshold_min   // channel: [mandatory] [ val(mitosalt_del_threshold_min) ]
        ch_deletion_threshold_max   // channel: [mandatory] [ val(mitosalt_del_threshold_max) ]
        ch_breakthreshold           // channel: [mandatory] [ val(mitosalt_break_threshold) ]
        ch_cluster_threshold        // channel: [mandatory] [ val(mitosalt_cluster_threshold) ]
        ch_breakspan                // channel: [mandatory] [ val(mitosalt_break_span) ]
        ch_sizelimit                // channel: [mandatory] [ val(mitosalt_size_limit) ]
        ch_hplimit                  // channel: [mandatory] [ val(mitosalt_hp_limit) ]
        ch_flank                    // channel: [mandatory] [ val(mitosalt_flank) ]
        ch_split_distance_threshold // channel: [mandatory] [ val(mitosalt_split_dist_threshold) ]
        ch_subdepth                 // channel: [mandatory] [ val(mitosalt_depth) ]
	ch_mito_name                // channel: [mandatory] [ val(mito_name) ]
	ch_exclude                  // channel: [mandatory] [ val(mitosalt_exclude) ]

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
		ch_exclude,                
	        ch_score_threshold,         
	        ch_evalue_threshold,        
	        ch_split_length,            
	        ch_paired_distance,         
	        ch_deletion_threshold_min,  
	        ch_deletion_threshold_max,  
	        ch_breakthreshold,          
	        ch_cluster_threshold,       
	        ch_breakspan,               
	        ch_sizelimit,               
	        ch_hplimit,                 
	        ch_flank,                   
	        ch_split_distance_threshold
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
