process CHANGE_NAME {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path( "*.${file_type}"), emit: file

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    file_type = task.ext.file_type ?: input_file.getExtension()
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mv \\
        $input_file \\
        ${prefix}.${file_type}
    """

    stub:
    def args = task.ext.args ?: ''
    file_type = task.ext.file_type ?: input_file.getExtension()
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${file_type}
    """
}
