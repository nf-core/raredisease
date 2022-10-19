//
// Merge and annotate MT
//

include { GATK4_MERGEVCFS as GATK4_MERGEVCFS_LIFT_UNLIFT_MT      } from '../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_VARIANTFILTRATION as GATK4_VARIANTFILTRATION_MT  } from '../../modules/nf-core/gatk4/variantfiltration/main'
include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS_MT                } from '../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX_TABIX_MT                          } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_NORM as REMOVE_DUPLICATES_MT                  } from '../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as TABIX_TABIX_MT2                         } from '../../modules/nf-core/tabix/tabix/main'
include { GATK4_MERGEVCFS as GATK4_MERGEVCFS_MT                  } from '../../modules/nf-core/gatk4/mergevcfs/main'
include { HMTNOTE as HMTNOTE_MT                                  } from '../../modules/nf-core/hmtnote/main'
include { TABIX_TABIX as TABIX_TABIX_MT3                         } from '../../modules/nf-core/tabix/tabix/main'
include { ENSEMBLVEP as ENSEMBLVEP_MT                            } from '../../modules/nf-core/ensemblvep/main'
include { HAPLOGREP2_CLASSIFY as HAPLOGREP2_CLASSIFY_MT          } from '../../modules/nf-core/haplogrep2/classify/main'

workflow MERGE_ANNOTATE_MT {
    take:
       vcf1          // channel: [ val(meta), path('*.vcf.gz') ]
       vcf2          // channel: [ val(meta), path('*.vcf.gz') ]
       fasta         // channel: [ genome.fasta ]
       dict          // channel: [ genome.dict ]
       fai           // channel: [ genome.fai ]
       vep_genome
       vep_cache_version
       vep_cache
       case_info      // channel: [ val(case_info) ]
       vep_genome
       vep_cache_version
       vep_cache

    main:
       ch_versions = Channel.empty()
       
        
        // Merging lifted and unlifted vcfs
        //vcfs=vcf1.buffer( size:1, skip:1 ).mix(vcf2.buffer( size:1, skip:1 ))
        // vcf1.map{ meta, vcf ->
        //     return [meta]}.set{ch_meta}
        // //ch_meta.view() 
        
        // vcf1.map{ meta, vcf ->
        //     return [vcf]}.set{ch_vcf1}
        // //ch_vcf1.view() 
        
        // vcf2.map{ meta, vcf ->
        //     return [vcf]}.set{ch_vcf2}
        // //ch_vcf2.view() 
        
        // ch_vcfs=ch_vcf1.concat(ch_vcf2).flatten().toList()
        //ch_vcfs=Channel.of(vcf1,vcf2).gtoupTuple().toList()
        //ch_vcfs.view()

        ch_merg_lift_unlif = vcf1.join(vcf2, by:[0])

        //ch_vcfs = ch_merg_lift_unlif.map{meta, vcf -> [[it[0]], [it[1].concat(it[2]).toList().flatten()]]}.groupTuple()
        //ch_vcfs = Channel.of(vcf1,vcf2).groupTuple().groupKey(meta)
        vcf1.view()
        ch_vcfs =Channel.empty().mix(vcf1,vcf2).map{meta, vcf -> [meta, groupKey(vcf,2)]}.groupTuple()
        //ch_vcfs =Channel.of(vcf1,vcf2).map{meta, vcf -> [meta, groupKey(vcf,2)]}.groupTuple()
        ch_vcfs.view()

        GATK4_MERGEVCFS_LIFT_UNLIFT_MT( ch_vcfs, dict)
        ch_versions = ch_versions.mix(GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.versions.first())
        
        // Filtering Mutect calls

        // ch_lift_unlift = GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.vcf.join(GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.tbi)
        // ch_to_filt = ch_lift_unlift.map {
        //     meta, vcf, tbi ->
        //         return [meta, vcf, tbi, [], [], [], [], []]}

        // GATK4_FILTERMUTECTCALLS_MT( ch_to_filt, 
        //     fasta, 
        //     fai, 
        //     dict )
        // ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS_MT.out.versions.first())
        
        // Filtering Variants
        ch_filt_vcf = GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.vcf.join(GATK4_MERGEVCFS_LIFT_UNLIFT_MT.out.tbi, by:[0])
        GATK4_VARIANTFILTRATION_MT(ch_filt_vcf, 
            fasta, 
            fai, 
            dict)
        ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION_MT.out.versions.first())
        
        // Spliting multiallelic calls
        ch_in_split=GATK4_VARIANTFILTRATION_MT.out.vcf.join( GATK4_VARIANTFILTRATION_MT.out.tbi, by:[0])
        SPLIT_MULTIALLELICS_MT (ch_in_split, fasta)
        ch_versions = ch_versions.mix(SPLIT_MULTIALLELICS_MT.out.versions)
        
        // Removing duplicates
        TABIX_TABIX_MT(SPLIT_MULTIALLELICS_MT.out.vcf)
        ch_in_remdup = SPLIT_MULTIALLELICS_MT.out.vcf.join(TABIX_TABIX_MT.out.tbi)
        REMOVE_DUPLICATES_MT(ch_in_remdup, fasta)
        ch_versions = ch_versions.mix(REMOVE_DUPLICATES_MT.out.versions)

        //TABIX_TABIX_MT2(REMOVE_DUPLICATES_MT.out.vcf)
        
        // Merging vcfs of different family members
        REMOVE_DUPLICATES_MT.out
            .vcf
            .collect{it[1]}
            .toList()
            .set { file_list }
        case_info
            .combine(file_list)
            .set { ch_mergvcf }
        GATK4_MERGEVCFS_MT( ch_mergvcf, dict)
        ch_versions = ch_versions.mix(GATK4_MERGEVCFS_MT.out.versions.first())
        
        
        // Annotating with Hmtnote
        //HMTNOTE_MT(GATK4_VARIANTFILTRATION_MT.out.vcf)
        //ch_versions = ch_versions.mix(HMTNOTE_MT.out.versions.first())

        //TABIX_TABIX_MT2(HMTNOTE_MT.out.vcf)
        //ch_versions = ch_versions.mix(TABIX_TABIX_MT2.out.versions.first())
        //ch_vep_in_mt=HMTNOTE_MT.out.vcf.join( TABIX_TABIX_MT2.out.tbi, by:[0])
        
        // Annotating with ensembl Vep
        ENSEMBLVEP_MT( GATK4_MERGEVCFS_MT.out.vcf,
            vep_genome,
            "homo_sapiens",
            vep_cache_version,
            vep_cache,
            fasta,
            [])
        ch_versions = ch_versions.mix(ENSEMBLVEP_MT.out.versions.first())
        
        // Running haplogrep2
        TABIX_TABIX_MT3(ENSEMBLVEP_MT.out.vcf_gz)
        HAPLOGREP2_CLASSIFY_MT(ENSEMBLVEP_MT.out.vcf_gz, "vcf.gz")
        ch_versions = ch_versions.mix(HAPLOGREP2_CLASSIFY_MT.out.versions.first())


    emit:
        haplog   = HAPLOGREP2_CLASSIFY_MT.out.txt
        vcf      = ENSEMBLVEP_MT.out.vcf_gz
        tbi      = TABIX_TABIX_MT3.out.tbi
        report   = ENSEMBLVEP_MT.out.report
        versions = ch_versions // channel: [ versions.yml ]
}