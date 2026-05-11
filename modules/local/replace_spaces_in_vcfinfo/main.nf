process REPLACE_SPACES_IN_VCFINFO {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2fdb3ffa4bac62e98ec062e991217df2ab96e67bbd45a502bd25ff2effecb96e/data':
        'community.wave.seqera.io/library/htslib_python:d1e4474cbf76f4e9' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*_reformatted.vcf"), emit: vcf
    tuple val("${task.process}"), val('replace_spaces_in_vcfinfo'), val("1.0"), topic: versions, emit: versions_replace_spaces_in_vcfinfo
    tuple val("${task.process}"), val('python'), eval("python -V | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    replace_spaces_in_vcfinfo.py --input ${input}
    """

    stub:
    """
    replace_spaces_in_vcfinfo.py --input ${input}
    """
}
