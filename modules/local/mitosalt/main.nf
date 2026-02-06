process MITOSALT {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dc2631dac526622d28ea1109b0f714d536606d0e5b3b85fe24407c8206e7e6b6/data':
        'community.wave.seqera.io/library/mitosalt:1.1.1--5fd87ac48a683358' }"


    input:
    tuple val(meta), path(reads)
    path msconfig
    path chrsizes
    tuple val(meta3), path(genomefai)
    tuple val(meta2), path(hisat2index)
    tuple val(meta5), path(mtfai)
    tuple val(meta6), path(mtfasta)
    tuple val(meta4), path(lastindex)

    output:
    tuple val(meta), path("*breakpoint") , emit: breakpoint
    tuple val(meta), path("*cluster")    , emit: cluster
    path "versions.yml"                  , emit: versions

    script:
    def VERSION = "1.1.1" // from perl script, unlikely to be updated
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $msconfig | sed "s/threads = 1/threads = ${task.cpus}/" > new-${msconfig}
    mkdir -p log indel bam tab bw plot
    mitosalt new-${msconfig} $reads $prefix
    mv indel/*.breakpoint ${prefix}.breakpoint
    mv indel/*.cluster ${prefix}.cluster

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitosalt: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = "1.1.1"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $msconfig | sed "s/threads = 1/threads = ${task.cpus}/" > new-${msconfig}
    touch ${prefix}.breakpoint
    touch ${prefix}.cluster

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitosalt: $VERSION
    END_VERSIONS
    """

}
