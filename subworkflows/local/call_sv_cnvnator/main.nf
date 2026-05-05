//
// A subworkflow to call CNVs using cnvnator
//

include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_CNVNATOR } from '../../../modules/nf-core/bcftools/view/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_CALL      } from '../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_HIST      } from '../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_PARTITION } from '../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_RD        } from '../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CNVNATOR as CNVNATOR_STAT      } from '../../../modules/nf-core/cnvnator/cnvnator/main.nf'
include { CNVNATOR_CONVERT2VCF                    } from '../../../modules/nf-core/cnvnator/convert2vcf/main.nf'
include { SPLIT_CHR                               } from '../../../modules/local/split_chr/main.nf'
include { SVDB_MERGE as SVDB_MERGE_CNVNATOR       } from '../../../modules/nf-core/svdb/merge/main'
include { TABIX_BGZIPTABIX as INDEX_CNVNATOR      } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow CALL_SV_CNVNATOR {
    take:
        ch_bam_bai   // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
        ch_fasta     // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_case_info // channel: [mandatory] [ val(case_info) ]

    main:

        SPLIT_CHR (ch_fasta)
        CNVNATOR_RD ( ch_bam_bai, [[:],[]], [[:],[]], [[:],[]], "rd" )
        CNVNATOR_HIST ( [[:],[],[]], CNVNATOR_RD.out.root, SPLIT_CHR.out.output, [[:],[]], "his" )
        CNVNATOR_STAT ( [[:],[],[]], CNVNATOR_HIST.out.root, [[:],[]], [[:],[]], "stat" )
        CNVNATOR_PARTITION ( [[:],[],[]], CNVNATOR_STAT.out.root, [[:],[]], [[:],[]], "partition" )
        CNVNATOR_CALL ( [[:],[],[]], CNVNATOR_PARTITION.out.root, [[:],[]], [[:],[]], "call" )
        CNVNATOR_CONVERT2VCF (CNVNATOR_CALL.out.tab)
        INDEX_CNVNATOR (CNVNATOR_CONVERT2VCF.out.vcf)
        BCFTOOLS_VIEW_CNVNATOR (INDEX_CNVNATOR.out.gz_index, [], [], []).vcf
            .collect{_meta, vcf -> vcf}
            .toList()
            .set { vcf_file_list }

        ch_case_info
            .combine(vcf_file_list)
            .set { merge_input_vcfs }

        SVDB_MERGE_CNVNATOR ( merge_input_vcfs, [], true )

    emit:
        vcf        = SVDB_MERGE_CNVNATOR.out.vcf  // channel: [ val(meta), path(*.tar.gz) ]
}
