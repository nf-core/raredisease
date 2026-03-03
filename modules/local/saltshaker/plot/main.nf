process SALTSHAKER_PLOT {
    tag "$meta.id"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e93d703b195dd27cd920cee46669d3f51043216c12fd05168c937e93adf170e8/data':
        'community.wave.seqera.io/library/pip_saltshaker:e08e38a6d45f8f32' }"

    input:
    tuple val(meta), path(classify)

    output:
    tuple val(meta), path("*saltshaker.png"), emit: plot
    path "versions.yml"                     , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    saltshaker plot \\
        --prefix $prefix \\
        --input-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        saltshaker_plot: \$(echo \$(saltshaker plot 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        saltshaker_plot: \$(echo \$(saltshaker plot 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

}
