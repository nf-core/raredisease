process SALTSHAKER_CALL {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e93d703b195dd27cd920cee46669d3f51043216c12fd05168c937e93adf170e8/data':
        'community.wave.seqera.io/library/pip_saltshaker:e08e38a6d45f8f32' }"

    input:
    tuple val(meta), path(breakpoint)
    tuple val(meta), path(cluster)
    tuple val(meta2), path(mtfasta)
    val flank
    val hplimit
    val mito_length
    val ori_h_start
    val ori_h_end
    val ori_l_start
    val ori_l_end

    output:
    tuple val(meta), path("*_call_metadata.tsv"), emit: call
    path "versions.yml"                         , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    saltshaker call \\
        --prefix $prefix \\
        --output-dir . \\
        -r $mtfasta \\
        -c $cluster \\
        -p $breakpoint \\
        -f $flank \\
        -H $hplimit \\
        -g $mito_length \\
        --ori-h-start $ori_h_start \\
        --ori-h-end $ori_h_end \\
        --ori-l-start $ori_l_start \\
        --ori-l-end $ori_l_end \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        saltshaker_call: \$(echo \$(saltshaker call 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_call_metadata.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        saltshaker_call: \$(echo \$(saltshaker call 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

}
