//
// A variant caller workflow for deepvariant
//

params.split_multiallelics_options = [:]
params.rm_duplicates_options = [:]
params.deepvariant_options = [:]
params.glnexus_options = [:]

def split_multiallelics_glnexus           = params.split_multiallelics_options.clone()
split_multiallelics_glnexus.publish_files = "false"

def rm_duplicates_glnexus                 = params.rm_duplicates_options.clone()
rm_duplicates_glnexus.publish_dir         = "glnexus/"
rm_duplicates_glnexus.suffix              = "_split_rmdup"


include { BCFTOOLS_NORM as SPLIT_MULTIALLELICS } from '../../modules/nf-core/modules/bcftools/norm/main'  addParams( options: split_multiallelics_glnexus )
include { BCFTOOLS_NORM as REMOVE_DUPLICATES } from '../../modules/nf-core/modules/bcftools/norm/main'  addParams( options: rm_duplicates_glnexus )
include { DEEPVARIANT } from '../../modules/local/deepvariant/main'  addParams( options: params.deepvariant_options )
include { GLNEXUS } from '../../modules/nf-core/modules/glnexus/main'  addParams( options: params.glnexus_options )

workflow DEEPVARIANT_CALLER {
    take:
        bam // channel: [ val(meta), path(bam), path(bai) ]
        fasta // path(fasta)
        fai // path(fai)
        sample // channel: [ sample, sex, phenotype, paternal_id, maternal_id, case_id ]

    main:
        DEEPVARIANT ( bam, fasta, fai )
        DEEPVARIANT.out.gvcf.collect{it[1]}
            .toList()
            .set { file_list }

        //retrieve case id for glnexus and store it in a new channel called case_meta
        sample
            .first()
            .map{
                it ->
                    new_sample_meta = it.clone()
                    new_sample_meta.id = new_sample_meta.case_id
                    [ [ 'id':new_sample_meta.id ] ] }
            .set {case_meta}

        //Combine case meta with the list of gvcfs
        case_meta.combine(file_list)
            .set { ch_gvcfs }
        GLNEXUS ( ch_gvcfs )
        SPLIT_MULTIALLELICS (GLNEXUS.out.bcf, fasta)
        REMOVE_DUPLICATES (SPLIT_MULTIALLELICS.out.vcf, fasta)

    emit:
        vcf                         = REMOVE_DUPLICATES.out.vcf

        // Collect versions
        deepvariant_version         = DEEPVARIANT.out.version
        glnexus_version             = GLNEXUS.out.version

}
