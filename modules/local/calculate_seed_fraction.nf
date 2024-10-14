process CALCULATE_SEED_FRACTION {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(cov)
    val rd
    val seed

    output:
    tuple val(meta), path("seedfrac.csv"), emit: csv
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MT_COVERAGE=`awk '{cov += \$3}END{ if (NR > 0) print cov / NR }' $cov`

    python -c "import os;print('%0.6f' % ($seed+ $rd/float(os.environ['MT_COVERAGE'])))" >seedfrac.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        calculate_seed_fraction: v1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch seedfrac.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        calculate_seed_fraction: v1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
