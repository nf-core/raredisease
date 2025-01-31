process GET_CHROM_SIZES {
    tag "$fai"
    label 'process_single'

    conda "conda-forge::coreutils=8.31"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    input:
    tuple val(meta), path(fai)

    output:
    path '*.sizes'     , emit: sizes
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = "8.31"
    """
    cut -f 1,2 $fai > ${fai}.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = "8.31"
    """
    touch ${fai}.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """

}
