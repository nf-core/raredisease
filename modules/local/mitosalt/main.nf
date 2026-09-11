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
    tuple val(meta2), path(genomefai)
    tuple val(meta3), path(mtfai)
    tuple val(meta4), path(mtfasta)
    tuple val(meta5), path(lastindex)

    output:
    tuple val(meta), path("*breakpoint") , emit: breakpoint
    tuple val(meta), path("*cluster")    , emit: cluster
    tuple val("${task.process}"), val('mitosalt'), val("1.1.1"), topic: versions, emit: versions_mitosalt

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    // Heap for the reformat.sh (BBMap) JVM step, sized to what this task was actually
    // given instead of the tool's own hardcoded -Xmx100g, which reserves more virtual
    // memory than any real allocation provides and fails JVM startup outright.
    def javamem = Math.max(1, (task.memory.toGiga() * 0.8) as int)
    """
    cat $msconfig | sed "s/threads = 1/threads = ${task.cpus}/" > new-${msconfig}
    echo "javamem = ${javamem}" >> new-${msconfig}
    mkdir -p log indel bam tab bw plot
    MitoSAlt1.1.1.pl new-${msconfig} $reads $prefix

    if grep -qiE "insufficient memory for the Java Runtime Environment|Native memory allocation \\(mmap\\) failed|Could not reserve enough space" log/${prefix}.log; then
        echo "ERROR: Fatal error detected in MitoSAlt log" >&2
        exit 137
    fi

    mv indel/*.breakpoint ${prefix}.breakpoint
    mv indel/*.cluster ${prefix}.cluster
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $msconfig | sed "s/threads = 1/threads = ${task.cpus}/" > new-${msconfig}
    touch ${prefix}.breakpoint
    echo 'cluster' > ${prefix}.cluster
    """

}
