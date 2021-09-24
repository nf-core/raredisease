//
// A variant caller workflow for deepvariant
//

params.deepvariant_options = [:]

include { DEEPVARIANT } from '../../modules/local/deepvariant/main'  addParams( options: params.deepvariant_options )
include { GLNEXUS } from '../../modules/nf-core/modules/glnexus/main'  addParams( options: params.deepvariant_options )

workflow DEEPVARIANT_CALLER {
    take:
        bam // channel: [ val(meta), path(bam), path(bai) ]
        fasta // path(fasta)
        fai // path(fai)
        sample // channel: [ sample, sex, phenotype, paternal_id, maternal_id, case_id ]

    main:
        DEEPVARIANT ( bam, fasta, fai )
        file_list = DEEPVARIANT.out.gvcf.collect{it[1]}.toList()

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
        case_meta.mix(file_list)
            .collect()
            .set { ch_gvcfs }
        GLNEXUS ( ch_gvcfs )

    emit:
        vcf                         = GLNEXUS.out.bcf

        // Collect versions
        deepvariant_version         = DEEPVARIANT.out.version
        glnexus_version             = GLNEXUS.out.version

}
