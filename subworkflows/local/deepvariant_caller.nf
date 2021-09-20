//
// A variant caller workflow for deepvariant
//

params.deepvariant_options = [:]

include { DEEPVARIANT } from '../../modules/local/deepvariant/main'  addParams( options: params.deepvariant_options )

workflow DEEPVARIANT_CALLER {
    take:
        bam // channel: [ val(meta), path(bam), path(bai) ]
        fasta // path(fasta)
        fai // path(fai)

    main:
        DEEPVARIANT ( bam, fasta, fai )

    emit:
        vcf                         = DEEPVARIANT.out.vcf
        gvcf                        = DEEPVARIANT.out.gvcf

        // Collect versions
        deepvariant_version         = DEEPVARIANT.out.version
}
