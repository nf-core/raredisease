//
// A variant caller workflow for deepvariant
//

params.split_multiallelics_options = [:]
params.rm_duplicates_options = [:]
params.deepvariant_options = [:]
params.glnexus_options = [:]
params.tabix_options = [:]

def split_multiallelics_glnexus           = params.split_multiallelics_options.clone()
split_multiallelics_glnexus.publish_files = "false"

def rm_duplicates_glnexus                 = params.rm_duplicates_options.clone()
rm_duplicates_glnexus.publish_dir         = "glnexus/"
rm_duplicates_glnexus.suffix              = "_split_rmdup"

def tabix_glnexus                         = params.tabix_options.clone()
tabix_glnexus.publish_dir                 = "glnexus/"

include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS } from '../../modules/nf-core/modules/bcftools/norm/main'  addParams( options: split_multiallelics_glnexus )
include { BCFTOOLS_NORM as REMOVE_DUPLICATES } from '../../modules/nf-core/modules/bcftools/norm/main'  addParams( options: rm_duplicates_glnexus )
include { DEEPVARIANT } from '../../modules/local/deepvariant/main'  addParams( options: params.deepvariant_options )
include { GLNEXUS } from '../../modules/nf-core/modules/glnexus/main'  addParams( options: params.glnexus_options )
include { TABIX_TABIX as TABIX } from '../../modules/nf-core/modules/tabix/tabix/main'  addParams( options: tabix_glnexus)

workflow DEEPVARIANT_CALLER {
    take:
    bam          // channel: [ val(meta), path(bam), path(bai) ]
    fasta        // path(fasta)
    fai          // path(fai)
    ch_case_info // channel: [ case_id ]

    main:
    DEEPVARIANT ( bam, fasta, fai )
    DEEPVARIANT.out.gvcf.collect{it[1]}
        .toList()
        .set { file_list }

    //Combine case meta with the list of gvcfs
    ch_case_info.combine(file_list)
        .set { ch_gvcfs }
    GLNEXUS ( ch_gvcfs )
    SPLIT_MULTIALLELICS (GLNEXUS.out.bcf, fasta)
    REMOVE_DUPLICATES (SPLIT_MULTIALLELICS.out.vcf, fasta)
    TABIX (REMOVE_DUPLICATES.out.vcf)

    emit:
    vcf                         = REMOVE_DUPLICATES.out.vcf
    tbi                         = TABIX.out.tbi

    // Collect versions
    deepvariant_version         = DEEPVARIANT.out.version
    glnexus_version             = GLNEXUS.out.version

}
