process MITOSALT {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mitosalt%3A1.1.1--hdfd78af_3':
        'community.wave.seqera.io/library/bbmap_bedtools_bioconductor-biostrings_bioconductor-pwalign_pruned:11434f3b6a01596d' }"


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
