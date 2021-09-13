//
// A variant caller workflow for deepvariant
//

params.deepvariant_options = [:]

include { DEEPVARIANT } from '../../modules/local/depvariant/main'  addParams( options: params.deepvariant_options )

workflow DEEPVARIANT_CALLER {
    take:
        bam // channel: [ val(meta), path(bam), path(bai) ]
        fasta // channel: [ path(fasta), path(fai) ]

    main:
        DEEPVARIANT ( bam, fasta )

    emit:
        vcf                         = DEEPVARIANT.out.vcf
        gvcf                        = DEEPVARIANT.out.gvf

        // Collect versions
        deepvariant_version         = DEEPVARIANT.out.version
}
