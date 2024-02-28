process RENAME_ALIGN_FILES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=8.31"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--0' :
        'biocontainers/gnu-wget:1.18--0' }"

    input:
    tuple val(meta), path(input)
    val(extension)

    output:
    path("*.{bam,bai}"), emit: output
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ln -s $input ${meta.sample}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ln: \$(echo \$(ln --version 2>&1 | head -n 1 | cut -d ' ' -f4))
    END_VERSIONS
    """
}
