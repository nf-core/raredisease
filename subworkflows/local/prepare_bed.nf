//
// Prepare reference bed files
//

include { GATK4_BEDTOINTERVALLIST as GATK_BILT } from '../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_INTERVALLISTTOOLS as GATK_ILT  } from '../../modules/nf-core/gatk4/intervallisttools/main'
include { CAT_CAT as CAT_CAT_BAIT              } from '../../modules/nf-core/cat/cat/main'
include { TABIX_TABIX as TABIX_PT              } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIPTABIX as TABIX_PBT        } from '../../modules/nf-core/tabix/bgziptabix/main'

workflow CHECK_BED {
    take:
        bed                // file: bed file
        seq_dictionary     // path: sequence_dictionary

    main:
        tab_out = Channel.empty()
        ch_versions = Channel.empty()

        if (bed) {
            bed_file = file(bed)
            id       = bed.split('/')[-1]
            ch_bed   = Channel.fromList([[['id':id], bed_file]])

            if ( bed.endsWith(".gz") && file(bed, checkIfExists:true) ) {
                tbi_out     = TABIX_PT (ch_bed).tbi
                tab_out     = ch_bed.join(tbi_out)
                ch_versions = ch_versions.mix(TABIX_PT.out.versions)
            } else if ( file(bed, checkIfExists:true) ) {
                tab_out     = TABIX_PBT (ch_bed).gz_tbi
                ch_versions = ch_versions.mix(TABIX_PBT.out.versions)
            }

            interval_list = GATK_BILT (ch_bed, seq_dictionary).interval_list
            ch_versions   = ch_versions.mix(GATK_BILT.out.versions)

            GATK_ILT(interval_list)
            ch_versions   = ch_versions.mix(GATK_ILT.out.versions)

            GATK_ILT.out
                .interval_list
                .collect{ it[1] }
                .set { ch_bait_intervals_split }

            ch_bait_intervals_split
                .map { it -> it[0]
                    .toString()
                    .split("_split")[0]
                    .split("/")[-1] + "_bait.intervals_list"
                }
                .flatten()
                .concat(ch_bait_intervals_split)
                .toList()
                .map {
                    id, bait ->
                        return [['id':id], bait]
                }
                .set { ch_bait_intervals_cat_in }

            CAT_CAT_BAIT ( ch_bait_intervals_cat_in )
                .file_out
                .map {
                    id, file ->
                        return [file]
                }
                .set { ch_bait_intervals_cat_out }

        }

    emit:
        bed              = tab_out
        target_intervals = interval_list.collect{it[1]}
        bait_intervals   = ch_bait_intervals_cat_out
        versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
