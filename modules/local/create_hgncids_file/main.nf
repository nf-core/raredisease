process CREATE_HGNCIDS_FILE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2fdb3ffa4bac62e98ec062e991217df2ab96e67bbd45a502bd25ff2effecb96e/data':
        'community.wave.seqera.io/library/htslib_python:d1e4474cbf76f4e9' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*_reformatted.txt"), emit: txt
    tuple val("${task.process}"), val('create_hgncids_file'), val("1.0"), topic: versions, emit: versions_create_hgncids_file
    tuple val("${task.process}"), val('python'), eval("python -V | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    create_hgncids_file.py --input ${input} --meta-id ${meta.id}

    """

    stub:
    """
    create_hgncids_file.py --input ${input} --meta-id ${meta.id}

    """
}
