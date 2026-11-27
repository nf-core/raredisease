//
// Merge and normalize MT variants
//

include { ADD_VARCALLER_TO_BED                                  } from '../../../modules/local/add_varcallername_to_bed'
include { BCFTOOLS_ANNOTATE                                     } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_MT                   } from '../../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES_MT                 } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_MT               } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_POSTMERGE_MT     } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_VIEW as KEEP_PASS_MT                         } from '../../../modules/nf-core/bcftools/view/main'
include { GATK4_MERGEVCFS as GATK4_MERGEVCFS_LIFT_UNLIFT_MT     } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_VARIANTFILTRATION as GATK4_VARIANTFILTRATION_MT } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { KEEP_PROBAND_VARIANTS                                 } from '../../../modules/local/keep_proband_variants/main'
include { PICARD_LIFTOVERVCF                                    } from '../../../modules/nf-core/picard/liftovervcf/main'
include { TABIX_TABIX as TABIX_ANNOTATE                         } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TABIX_MERGE                      } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TABIX_MT                         } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TABIX_MT2                        } from '../../../modules/nf-core/tabix/tabix/main'

workflow POSTPROCESS_MT_CALLS {
    take:
        ch_case_info         // channel: [mandatory] [ val(case_info) ]
        ch_foundin_header    // channel: [mandatory] [ path(header) ]
        ch_genome_chrsizes   // channel: [mandatory] [ path(chrsizes) ]
        ch_genome_dictionary // channel: [mandatory] [ val(meta), path(dict) ]
        ch_genome_fai        // channel: [mandatory] [ val(meta), path(fai) ]
        ch_genome_fasta      // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_mt_vcf            // channel: [mandatory] [ val(meta), path(vcf) ]
        ch_mtshift_backchain // channel: [mandatory] [ val(meta), path(backchain) ]
        ch_mtshift_vcf       // channel: [mandatory] [ val(meta), path(vcf) ]

    main:
        // LIFTOVER SHIFTED VCF TO REFERENCE MT POSITIONS
        PICARD_LIFTOVERVCF (
            ch_mtshift_vcf,
            ch_genome_dictionary,
            ch_genome_fasta,
            ch_mtshift_backchain,
        )

        ch_vcfs = ch_mt_vcf
            .join(PICARD_LIFTOVERVCF.out.vcf_lifted, remainder: true)
            .map{ meta, vcf1, vcf2 ->
                [meta, [vcf1, vcf2]]
            }

        GATK4_MERGEVCFS_LIFT_UNLIFT_MT( ch_vcfs, ch_genome_dictionary)

        // Filtering Variants
        ch_filt_vcf = GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.vcf
            .join(GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        GATK4_VARIANTFILTRATION_MT (ch_filt_vcf, ch_genome_fasta, ch_genome_fai, ch_genome_dictionary, [[:],[]])

        // Keep only PASS records in each sample before merging, so that a filtered call in one
        // sample can neither remove nor rescue a site in another sample once the VCFs are merged
        ch_in_pass = GATK4_VARIANTFILTRATION_MT.out.vcf
            .join(GATK4_VARIANTFILTRATION_MT.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        KEEP_PASS_MT (ch_in_pass, [], [], [])

        // Spliting multiallelic calls
        ch_in_split = KEEP_PASS_MT.out.vcf
            .join(KEEP_PASS_MT.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        SPLIT_MULTIALLELICS_MT (ch_in_split, ch_genome_fasta)
        TABIX_TABIX_MT(SPLIT_MULTIALLELICS_MT.out.vcf)

        // Removing duplicates and merging if there is more than one sample
        ch_in_remdup = SPLIT_MULTIALLELICS_MT.out.vcf
            .join(TABIX_TABIX_MT.out.index, failOnMismatch:true, failOnDuplicate:true)
        REMOVE_DUPLICATES_MT(ch_in_remdup, ch_genome_fasta)
        TABIX_TABIX_MT2(REMOVE_DUPLICATES_MT.out.vcf)

        file_list_vcf = REMOVE_DUPLICATES_MT.out.vcf
            .map{ _meta, vcf -> vcf}
            .toSortedList{a, b -> a.name <=> b.name}
            .toList()

        file_list_tbi = TABIX_TABIX_MT2.out.index
            .map{ _meta, vcf -> vcf}
            .toSortedList{a, b -> a.name <=> b.name}
            .toList()

        ch_rem_dup_vcf_tbi = ch_case_info
            .combine(file_list_vcf)
            .combine(file_list_tbi)

        ch_case_vcf = ch_rem_dup_vcf_tbi.branch {
            meta, vcf, tbi ->
                single: vcf.size() == 1
                    return [meta, vcf]
                multiple: vcf.size() > 1
                    return [meta, vcf, tbi]
            }

        BCFTOOLS_MERGE_MT(
            ch_case_vcf.multiple.map { meta, vcf, tbi -> return [meta, vcf, tbi, []] },
            ch_genome_fasta.join(ch_genome_fai, failOnMismatch:true, failOnDuplicate:true).collect()
        )

        SPLIT_MULTIALLELICS_POSTMERGE_MT(
            BCFTOOLS_MERGE_MT.out.vcf.join(BCFTOOLS_MERGE_MT.out.index, failOnMismatch:true, failOnDuplicate:true).collect(),
            ch_genome_fasta
        )

        // Keep only sites carried by an affected sample: mtDNA is maternally inherited, so a variant
        // seen only in an unaffected relative should not become a candidate in the family VCF.
        // Only relevant for merged multi-sample cases, and only possible when the case has affected samples.
        ch_merged_vcf = SPLIT_MULTIALLELICS_POSTMERGE_MT.out.vcf
            .branch { meta, _vcf ->
                with_probands: meta.probands
                without_probands: true
            }

        KEEP_PROBAND_VARIANTS(ch_merged_vcf.with_probands)

        ch_addfoundintag_in = KEEP_PROBAND_VARIANTS.out.vcf
            .mix(ch_merged_vcf.without_probands)
            .mix(ch_case_vcf.single)

        TABIX_TABIX_MERGE(ch_addfoundintag_in)

        ch_varcallerinfo = ch_genome_chrsizes.flatten().map{chromsizes ->
            return [[id:'mutect2'], chromsizes]
            }

        ADD_VARCALLER_TO_BED (ch_varcallerinfo)
        ch_varcallerbed = ADD_VARCALLER_TO_BED.out.gz_tbi
            .map{_meta,bed,tbi -> return [bed, tbi]}

        ch_annotate_in = ch_addfoundintag_in
            .join(TABIX_TABIX_MERGE.out.index)
            .combine(ch_varcallerbed)
            .combine(ch_foundin_header)
            .map { meta, vcf, vcf_tbi, bed, bed_tbi, hdr -> return [meta, vcf, vcf_tbi, bed, bed_tbi, [], hdr, []] }

        BCFTOOLS_ANNOTATE(ch_annotate_in)

        TABIX_ANNOTATE(BCFTOOLS_ANNOTATE.out.vcf)

    emit:
        tbi = TABIX_ANNOTATE.out.index  // channel: [ val(meta), path(tbi) ]
        vcf = BCFTOOLS_ANNOTATE.out.vcf // channel: [ val(meta), path(vcf) ]
}
