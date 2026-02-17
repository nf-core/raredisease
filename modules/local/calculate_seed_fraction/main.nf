process CALCULATE_SEED_FRACTION {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2fdb3ffa4bac62e98ec062e991217df2ab96e67bbd45a502bd25ff2effecb96e/data':
        'community.wave.seqera.io/library/htslib_python:d1e4474cbf76f4e9' }"

    input:
    tuple val(meta), path(cov)
    val rd
    val seed

    output:
    tuple val(meta), path("seedfrac.csv"), emit: csv
    tuple val("${task.process}"), val('calculate_seed_fraction'), val("1.0"), topic: versions, emit: versions_tabix
    tuple val("${task.process}"), val('python'), eval("python -V | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export MT_COVERAGE=`awk '{cov += \$3}END{ if (NR > 0) print cov / NR }' $cov`

    python -c "import os;print('%0.6f' % ($seed+ $rd/float(os.environ['MT_COVERAGE'])))" >seedfrac.csv
    """

    stub:
    """
    touch seedfrac.csv
    """
}
