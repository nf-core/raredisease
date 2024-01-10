//
// A variant caller workflow for GATK's GermlinceCNVCaller
//

include { GATK4_COLLECTREADCOUNTS             } from '../../../modules/nf-core/gatk4/collectreadcounts/main.nf'
include { GATK4_DETERMINEGERMLINECONTIGPLOIDY } from '../../../modules/nf-core/gatk4/determinegermlinecontigploidy/main.nf'
include { GATK4_GERMLINECNVCALLER             } from '../../../modules/nf-core/gatk4/germlinecnvcaller/main.nf'
include { GATK4_POSTPROCESSGERMLINECNVCALLS   } from '../../../modules/nf-core/gatk4/postprocessgermlinecnvcalls/main.nf'
include { BCFTOOLS_VIEW                       } from '../../../modules/nf-core/bcftools/view/main'
include { TABIX_TABIX                         } from '../../../modules/nf-core/tabix/tabix/main'

workflow CALL_SV_GERMLINECNVCALLER {
    take:
        ch_bam_bai             // channel: [mandatory][ val(meta), path(bam), path(bai) ]
        ch_fasta               // channel: [mandatory][ val(meta), path(ch_fasta_no_meta) ]
        ch_fai                 // channel: [mandatory][ val(meta), path(ch_fai) ]
        ch_readcount_intervals // channel: [mandatory][ path(intervals) ]
        ch_genome_dictionary   // channel: [mandatory][ val(meta), path(ch_dict) ]
        ch_ploidy_model        // channel: [mandatory][ path(ch_ploidy_model) ]
        ch_gcnvcaller_model    // channel: [mandatory][ path(ch_gcnvcaller_model) ]

    main:
        ch_versions = Channel.empty()

        input = ch_bam_bai.combine( ch_readcount_intervals )

        GATK4_COLLECTREADCOUNTS ( input, ch_fasta, ch_fai, ch_genome_dictionary )

        GATK4_COLLECTREADCOUNTS.out.tsv
                .map({ meta, tsv -> return [meta, tsv, [], [] ]})
                .set{ch_dgcp_in}

        GATK4_DETERMINEGERMLINECONTIGPLOIDY ( ch_dgcp_in, ch_ploidy_model, [] )

        GATK4_COLLECTREADCOUNTS.out.tsv
                .join(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls)
                .combine(ch_gcnvcaller_model)
                .map({ meta, tsv, calls, meta2, model -> return [meta, tsv, [], calls, model ]})
                .set{ch_gcnvc_in}

        GATK4_GERMLINECNVCALLER ( ch_gcnvc_in )

        GATK4_GERMLINECNVCALLER.out.calls.toList()
            .flatMap {reduce_input(it)}
            .buffer (size: 2)
            .combine(ch_gcnvcaller_model.collect{it[1]}.toList())
            .join(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.calls)
            .set {ch_postproc_in}

        GATK4_POSTPROCESSGERMLINECNVCALLS ( ch_postproc_in )

        TABIX_TABIX(GATK4_POSTPROCESSGERMLINECNVCALLS.out.segments)
        GATK4_POSTPROCESSGERMLINECNVCALLS.out.segments
            .join(TABIX_TABIX.out.tbi, failOnMismatch:true)
            .set {ch_segments_in}
        // Filter out reference only (0/0) segments
        BCFTOOLS_VIEW (ch_segments_in , [], [], [] )

        ch_versions = ch_versions.mix(GATK4_COLLECTREADCOUNTS.out.versions)
        ch_versions = ch_versions.mix(GATK4_DETERMINEGERMLINECONTIGPLOIDY.out.versions)
        ch_versions = ch_versions.mix(GATK4_GERMLINECNVCALLER.out.versions)
        ch_versions = ch_versions.mix(GATK4_POSTPROCESSGERMLINECNVCALLS.out.versions)
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)
        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

    emit:
        genotyped_intervals_vcf          = GATK4_POSTPROCESSGERMLINECNVCALLS.out.intervals  // channel: [ val(meta), path(*.vcf.gz) ]
        genotyped_segments_vcf           = GATK4_POSTPROCESSGERMLINECNVCALLS.out.segments   // channel: [ val(meta), path(*.vcf.gz) ]
        genotyped_filtered_segments_vcf  = BCFTOOLS_VIEW.out.vcf                            // channel: [ val(meta), path(*.vcf.gz) ]
        denoised_vcf                     = GATK4_POSTPROCESSGERMLINECNVCALLS.out.denoised   // channel: [ val(meta), path(*.vcf.gz) ]
        versions                         = ch_versions                                      // channel: [ versions.yml ]
}

// This function groups calls with same meta for postprocessing.
def reduce_input (List gcnvoutput) {
    def dictionary  = [:]
    def reducedList = []
    for (int i = 0; i<gcnvoutput.size(); i++) {
            meta  = gcnvoutput[i][0]
            model = gcnvoutput[i][1]
        if(dictionary.containsKey(meta)) {
            dictionary[meta] += [model]
        } else {
            dictionary[meta]  = [model]
        }
    }

    for (i in dictionary) {
        reducedList.add(i.key)
        reducedList.add(i.value)
        }
    return reducedList
}
