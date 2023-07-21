process GET_CHROM_SIZES {
    tag "$fai"
    label 'process_single'

    conda "conda-forge::coreutils=8.31"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--0' :
        'biocontainers/gnu-wget:1.18--0' }"

    input:
    tuple val(meta), path(fai)

    output:
    path '*.sizes'     , emit: sizes
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cut -f 1,2 $fai > ${fai}.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cut: \$(echo \$(cut --help 2>&1 | head -n 1 | cut -f1,2 -d' '))
    END_VERSIONS
    """
}
