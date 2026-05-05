//
// Call SNV MT
//

include { GATK4_FILTERMUTECTCALLS as  GATK4_FILTERMUTECTCALLS_MT } from '../../../modules/nf-core/gatk4/filtermutectcalls/main'
include { GATK4_MUTECT2 as GATK4_MUTECT2_MT                      } from '../../../modules/nf-core/gatk4/mutect2/main'

workflow CALL_SNV_MT {
    take:
        ch_bam_bai    // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_dict       // channel: [mandatory] [ val(meta), path(dict) ]
        ch_fai        // channel: [mandatory] [ val(meta), path(fai) ]
        ch_fasta      // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_intervals  // channel: [mandatory] [ path(interval_list) ]

    main:

    ch_bam_bai_int = ch_bam_bai.combine(ch_intervals)

        GATK4_MUTECT2_MT (ch_bam_bai_int, ch_fasta, ch_fai.map{meta, fai -> return [meta,fai,[]]}, ch_dict, [], [], [], [], [],[])

        // Filter Mutect2 calls
        ch_mutect_vcf = GATK4_MUTECT2_MT.out.vcf.join(GATK4_MUTECT2_MT.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        ch_mutect_out = ch_mutect_vcf.join(GATK4_MUTECT2_MT.out.stats, failOnMismatch:true, failOnDuplicate:true)
        ch_to_filt    = ch_mutect_out.map {
                            meta, vcf, tbi, stats ->
                            return [meta, vcf, tbi, stats, [], [], [], []]
                        }

        GATK4_FILTERMUTECTCALLS_MT (ch_to_filt, ch_fasta, ch_fai, ch_dict)

    emit:
        filt_stats     = GATK4_FILTERMUTECTCALLS_MT.out.stats // channel: [ val(meta), path(tsv) ]
        stats          = GATK4_MUTECT2_MT.out.stats           // channel: [ val(meta), path(stats) ]
        tbi            = GATK4_FILTERMUTECTCALLS_MT.out.tbi   // channel: [ val(meta), path(tbi) ]
        vcf            = GATK4_FILTERMUTECTCALLS_MT.out.vcf   // channel: [ val(meta), path(vcf) ]
}
