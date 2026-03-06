//
// A variant caller workflow for deepvariant
//

include { BCFTOOLS_ANNOTATE                          } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_GL    } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES_GL      } from '../../../modules/nf-core/bcftools/norm/main'
include { DEEPVARIANT_RUNDEEPVARIANT as DEEPVARIANT  } from '../../../modules/nf-core/deepvariant/rundeepvariant/main'
include { GLNEXUS                                    } from '../../../modules/nf-core/glnexus/main'
include { TABIX_BGZIP                                } from '../../../modules/nf-core/tabix/bgzip/main'
include { ADD_VARCALLER_TO_BED                       } from '../../../modules/local/add_varcallername_to_bed'

workflow CALL_SNV_DEEPVARIANT {
    take:
        ch_bam_bai         // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_case_info       // channel: [mandatory] [ val(case_info) ]
        ch_foundin_header  // channel: [mandatory] [ path(header) ]
        ch_genome_chrsizes // channel: [mandatory] [ path(chrsizes) ]
        ch_genome_fai      // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_fasta    // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_par_bed         // channel: [optional] [ val(meta), path(bed) ]
        ch_target_bed      // channel: [mandatory] [ val(meta), path(bed), path(index) ]
        val_analysis_type  // boolean

    main:

        if (val_analysis_type.equals("wes")) {
            TABIX_BGZIP(ch_target_bed.map{meta, gzbed, _index -> return [meta, gzbed]})
            ch_bam_bai
                .combine (TABIX_BGZIP.out.output.map {_meta, bed -> return bed})
                .set { ch_deepvar_in }
        } else if (val_analysis_type.equals("wgs")) {
            ch_bam_bai
                .map { meta, bam, bai ->
                        return [meta, bam, bai, []] }
                .set { ch_deepvar_in }
        }

        DEEPVARIANT ( ch_deepvar_in, ch_genome_fasta, ch_genome_fai, [[],[]], ch_par_bed )
        DEEPVARIANT.out.gvcf
            .map{ _meta, gvcf -> gvcf}
            .toSortedList{a, b -> a.name <=> b.name}
            .toList()
            .set { ch_file_list }

        ch_case_info
            .combine(ch_file_list)
            .map {meta, gvcf -> return [meta, gvcf, []]}
            .set { ch_gvcfs }

        GLNEXUS ( ch_gvcfs, [[:],[]] )

        ch_split_multi_in = GLNEXUS.out.bcf
                            .map{ meta, bcf ->
                                    return [meta, bcf, []] }
        SPLIT_MULTIALLELICS_GL (ch_split_multi_in, ch_genome_fasta)

        ch_remove_dup_in = SPLIT_MULTIALLELICS_GL.out.vcf
                            .map{ meta, vcf ->
                                    return [meta, vcf, []] }
        REMOVE_DUPLICATES_GL (ch_remove_dup_in, ch_genome_fasta)

        ch_genome_chrsizes.flatten().map{chromsizes ->
            return [[id:'deepvariant'], chromsizes]
            }
            .set { ch_varcallerinfo }

        ADD_VARCALLER_TO_BED (ch_varcallerinfo).gz_tbi
            .map{_meta,bed,tbi -> return [bed, tbi]}
            .set{ch_varcallerbed}

        REMOVE_DUPLICATES_GL.out.vcf
            .join(REMOVE_DUPLICATES_GL.out.tbi)
            .combine(ch_varcallerbed)
            .map { it -> it + [[]] }
            .combine(ch_foundin_header)
            .map { it -> it + [[]] }
            .set { ch_annotate_in }

        BCFTOOLS_ANNOTATE(ch_annotate_in)

    emit:
        gvcf       = DEEPVARIANT.out.gvcf       // channel: [ val(meta), path(gvcf)]
        gvcf_tabix = DEEPVARIANT.out.gvcf_tbi   // channel: [ val(meta), path(gvcf_tbi)]
        tabix      = BCFTOOLS_ANNOTATE.out.tbi  // channel: [ val(meta), path(tbi) ]
        vcf        = BCFTOOLS_ANNOTATE.out.vcf  // channel: [ val(meta), path(vcf) ]
}
