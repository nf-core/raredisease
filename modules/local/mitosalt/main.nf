process MITOSALT {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/13/130779a0dd5a8d86441f262a16a5fb1dfe562125edf93646b53c893d982b5519/data':
        'community.wave.seqera.io/library/bbmap_bedtools_bioconductor-biostrings_bioconductor-pwalign_pruned:856c05081cbd8239' }"

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
    tuple val("${task.process}"), val('mitosalt'), val("1.1.1"), topic: versions, emit: versions_mitosalt

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $msconfig | sed "s/threads = 1/threads = ${task.cpus}/" > new-${msconfig}
    mkdir -p log indel bam tab bw plot
    MitoSAlt1.1.1.pl new-${msconfig} $reads $prefix
    mv indel/*.breakpoint ${prefix}.breakpoint
    mv indel/*.cluster ${prefix}.cluster
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $msconfig | sed "s/threads = 1/threads = ${task.cpus}/" > new-${msconfig}
    touch ${prefix}.breakpoint
    touch ${prefix}.cluster
    """

}
