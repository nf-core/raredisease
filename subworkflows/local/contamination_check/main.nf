//
// Subworkflow: Check sample contamination using GATK
//

include { GATK4_GETPILEUPSUMMARIES     } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION } from '../../../modules/nf-core/gatk4/calculatecontamination/main'

workflow CONTAMINATION_CHECK {
    
    take:
    ch_bam               // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta             // channel: [ val(meta), path(fasta) ]
    ch_fai               // channel: [ val(meta), path(fai) ]
    ch_dict              // channel: [ val(meta), path(dict) ]
    ch_contamination_vcf // channel: [ path(vcf), path(tbi) ]
    ch_intervals         // channel: [ path(bed) ] - only used for WES, empty for WGS

    main:
    ch_versions = Channel.empty()
    
    // Prepare BAM with intervals - conditionally based on analysis type
    // For WGS: intervals will be empty [], for WES: intervals will contain the BED file
    ch_bam_with_intervals = ch_bam
        .combine(ch_intervals.ifEmpty([[]]))
        .map { meta, bam, bai, bed -> 
            // If bed is an empty list, pass empty list; otherwise pass bed file
            def intervals = (bed instanceof List && bed.isEmpty()) ? [] : (bed ?: [])
            [ meta, bam, bai, intervals ]
        }
    
    // Separate VCF and TBI - collect them to make value channels
    ch_variants_vcf = ch_contamination_vcf
        .map { vcf, tbi -> vcf }
        .collect()
    
    ch_variants_tbi = ch_contamination_vcf
        .map { vcf, tbi -> tbi }
        .collect()
    
    // Run GetPileupSummaries
    GATK4_GETPILEUPSUMMARIES (
        ch_bam_with_intervals,  // [meta, bam, bai, intervals]
        ch_fasta,               // [meta2, fasta]
        ch_fai,                 // [meta3, fai]
        ch_dict,                // [meta4, dict]
        ch_variants_vcf,        // path(vcf)
        ch_variants_tbi         // path(tbi)
    )
    ch_versions = ch_versions.mix(GATK4_GETPILEUPSUMMARIES.out.versions.first())
    
    // Run CalculateContamination (tumor-only, no matched normal)
    // Format: [meta, pileup_table, matched_normal_table]
    // matched_normal_table is empty [] for tumor-only mode
    ch_contamination_input = GATK4_GETPILEUPSUMMARIES.out.table
        .map { meta, table -> [ meta, table, [] ] }
    
    GATK4_CALCULATECONTAMINATION (
        ch_contamination_input
    )
    ch_versions = ch_versions.mix(GATK4_CALCULATECONTAMINATION.out.versions.first())

    emit:
    contamination_table = GATK4_CALCULATECONTAMINATION.out.contamination
    segmentation_table  = GATK4_CALCULATECONTAMINATION.out.segmentation
    pileup_table        = GATK4_GETPILEUPSUMMARIES.out.table
    versions            = ch_versions
}