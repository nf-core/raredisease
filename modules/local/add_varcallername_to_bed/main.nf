process ADD_VARCALLER_TO_BED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(chromsizes)

    output:
    tuple val(meta), path("*.gz"), path("*.tbi"), emit: gz_tbi
    path  "versions.yml" ,                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def variant_caller = "${meta.id}"
    """
    awk '{print \$1"\t0\t"\$2"\t$variant_caller\"}' $chromsizes > ${variant_caller}.bed
    bgzip --threads ${task.cpus} -c $args ${variant_caller}.bed > ${prefix}.bed.gz
    tabix $args2 ${prefix}.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed.gz
    touch ${prefix}.bed.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
