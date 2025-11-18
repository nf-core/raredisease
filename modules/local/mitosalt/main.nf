#!/usr/bin/env nextflow

process MITOSALT {
    tag "MITOSALT"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/03/0311c283e73736be01c2cbd1ca93ae826c209d9733ffa6d2d4d2caa31e7464cc/data':
        'community.wave.seqera.io/library/bbmap_bedtools_bioconductor-biostrings_bioconductor-pwalign_pruned:11434f3b6a01596d' }"

    input:
    tuple val(meta), path(reads)
    path config
    path genome

    output:
    tuple val(meta), path("*breakpoint") , emit: breakpoint
    tuple val(meta), path("*cluster")    , emit: cluster
    path "versions.yml"                  , emit: versions

    script:
    def VERSION = "1.1.1" // from perl script, unlikely to be updated
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    MitoSAlt1.1.1.pl $config $reads $prefix
    mv indel/*.breakpoint mitosalt_${prefix}.breakpoint
    mv indel/*.cluster mitosalt_${prefix}.cluster

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitosalt: $VERSION
    END_VERSIONS
    """
}
