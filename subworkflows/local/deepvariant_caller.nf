//
// A variant caller workflow for deepvariant
//
include { getCaseMeta } from './functions'

params.deepvariant_options = [:]

include { DEEPVARIANT } from '../../modules/local/deepvariant/main'  addParams( options: params.deepvariant_options )
include { GLNEXUS } from '../../modules/nf-core/modules/glnexus/main'  addParams( options: params.deepvariant_options )

workflow DEEPVARIANT_CALLER {
    take:
        bam     // channel: [ val(meta), path(bam), path(bai) ]
        fasta   // path(fasta)
        fai     // path(fai)
        sample  // channel: [ sample, sex, phenotype, paternal_id, maternal_id, case_id ]

    main:
        //
        // Run DeepVariant and create a channel containing a list of gvcf files
        //
        DEEPVARIANT ( bam, fasta, fai )
        ch_dv_gvcfs = DEEPVARIANT.out.gvcf.collect{it[1]}.toList()

        //
        // Retrieve case id for glnexus and store it in a new channel called case_meta
        //
        ch_case_meta = getCaseMeta(sample)

        //
        //Combine case meta with the list of gvcfs
        //
        ch_gvcfs = ch_case_meta.combine(ch_dv_gvcfs)

        //
        // Run glnexus and normalize the output
        //
        GLNEXUS ( ch_gvcfs )

    emit:
        vcf                         = GLNEXUS.out.bcf

        // Collect versions
        deepvariant_version         = DEEPVARIANT.out.version
        glnexus_version             = GLNEXUS.out.version

}
