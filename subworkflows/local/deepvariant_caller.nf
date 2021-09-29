//
// A variant caller workflow for deepvariant
//

params.deepvariant_options = [:]
params.glnexus_options = [:]

include { DEEPVARIANT } from '../../modules/local/deepvariant/main'  addParams( options: params.deepvariant_options )
include { GLNEXUS } from '../../modules/nf-core/modules/glnexus/main'  addParams( options: params.glnexus_options )

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

    emit:
        vcf                         = GLNEXUS.out.bcf

        // Collect versions
        deepvariant_version         = DEEPVARIANT.out.version
        glnexus_version             = GLNEXUS.out.version

}
