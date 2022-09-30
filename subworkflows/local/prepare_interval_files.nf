//
// Prepare the interval list files for Picard's CollectWgsMetrics
//

include { GATK4_BEDTOINTERVALLIST as GATK_BILT_WG } from '../../modules/nf-core/modules/gatk4/bedtointervallist/main'
include { TABIX_TABIX as TABIX_PT_WG              } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_BGZIPTABIX as TABIX_PBT_WG        } from '../../modules/nf-core/modules/tabix/bgziptabix/main'
include { GATK4_BEDTOINTERVALLIST as GATK_BILT_Y  } from '../../modules/nf-core/modules/gatk4/bedtointervallist/main'
include { TABIX_TABIX as TABIX_PT_Y               } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_BGZIPTABIX as TABIX_PBT_Y         } from '../../modules/nf-core/modules/tabix/bgziptabix/main'

workflow PREPARE_INTERVAL {
    take:
        aligner        // [mandatory] params.aligner
        bed_wg         // bed file with autosomes, sex chromosomes, and mt
        bed_y          // bed file with Y chromosome
        seq_dictionary // sequence dictionary

    main:
        ch_tab_wg_out  = Channel.empty()
        ch_tab_y_out   = Channel.empty()
        ch_versions = Channel.empty()

        if ( aligner != "sentieon" ) { //these files are not necessary if one runs sentieon
            //Prepare the WG bed file
            bed_wg_file = file(bed_wg)
            id_wg       = bed_wg.split('/')[-1]
            ch_bed_wg   = Channel.fromList([[['id':id_wg], bed_wg_file]])

            if ( bed_wg.endsWith(".gz") && file(bed_wg, checkIfExists:true) ) {
                tbi_wg_out  = TABIX_PT_WG (ch_bed_wg).tbi
                ch_tab_wg_out  = ch_bed_wg.join(tbi_wg_out)
                ch_versions = ch_versions.mix(TABIX_PT_WG.out.versions)
            } else if ( file(bed_wg, checkIfExists:true) ) {
                ch_tab_wg_out  = TABIX_PBT_WG (ch_bed_wg).gz_tbi
                ch_versions = ch_versions.mix(TABIX_PBT_WG.out.versions)
            }

            ch_interval_list_wg = GATK_BILT_WG (ch_bed_wg, seq_dictionary).interval_list
            ch_versions      = ch_versions.mix(GATK_BILT_WG.out.versions)

            //Prepare the Y chr bed file
            bed_y_file = file(bed_y)
            id_y       = bed_y.split('/')[-1]
            ch_bed_y   = Channel.fromList([[['id':id_y], bed_y_file]])

            if ( bed_y.endsWith(".gz") && file(bed_y, checkIfExists:true) ) {
                tbi_y_out   = TABIX_PT_Y (ch_bed_y).tbi
                ch_tab_y_out   = ch_bed_y.join(tbi_y_out)
                ch_versions = ch_versions.mix(TABIX_PT_Y.out.versions)
            } else if ( file(bed_y, checkIfExists:true) ) {
                ch_tab_y_out   = TABIX_PBT_Y (ch_bed_y).gz_tbi
                ch_versions = ch_versions.mix(TABIX_PBT_Y.out.versions)
            }

            ch_interval_list_y = GATK_BILT_Y (ch_bed_y, seq_dictionary).interval_list
            ch_versions     = ch_versions.mix(GATK_BILT_Y.out.versions)
        }

    emit:
        bed_wg_gz    = ch_tab_wg_out
        intervals_wg = ch_interval_list_wg
        bed_y_gz     = ch_tab_y_out
        intervals_y  = ch_interval_list_y
        versions     = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}

