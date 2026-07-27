process SANITY_CHECK_VCFANNO_DATABASES {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2fdb3ffa4bac62e98ec062e991217df2ab96e67bbd45a502bd25ff2effecb96e/data':
        'community.wave.seqera.io/library/htslib_python:d1e4474cbf76f4e9' }"

    input:
    path(toml)
    path(databases)

    output:
    path("*_filtered.toml")                                                                                  , emit: toml
    tuple val("${task.process}"), val('sanity_check_vcfanno_databases'), val("1.0"), topic: versions         , emit: versions_sanity_check_vcfanno_databases
    tuple val("${task.process}"), val('python'), eval("python -V | sed 's/Python //'"), topic: versions      , emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    sanity_check_vcfanno_databases.py \\
        --toml ${toml} \\
        --databases ${databases}

    """

    stub:
    """
    touch ${toml.baseName}_filtered.toml

    """
}
