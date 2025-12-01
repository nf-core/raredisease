#!/usr/bin/env nextflow

process MITOSALT {
    tag "$meta.id"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/03/0311c283e73736be01c2cbd1ca93ae826c209d9733ffa6d2d4d2caa31e7464cc/data':
        'community.wave.seqera.io/library/bbmap_bedtools_bioconductor-biostrings_bioconductor-pwalign_pruned:11434f3b6a01596d' }"

//    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'https://depot.galaxyproject.org/singularity/mitosalt:1.1.1--hdfd78af_2':
//        'quay.io/biocontainers/mitosalt:1.1.1--hdfd78af_2' }"
    
//    container 'docker://ieduba/mitosalt:latest'

    input:
    tuple val(meta), path(reads)
    path config
    tuple val(meta2), path(hisat2index)
    tuple val(meta3), path(genomefai)
    tuple val(meta4), path(lastindex)
    tuple val(meta5), path(mtfai)
    path chrsizes
    tuple val(meta6), path(mtfasta)
    val mito_name

    output:
    tuple val(meta), path("*breakpoint") , emit: breakpoint
    tuple val(meta), path("*cluster")    , emit: cluster
    path "versions.yml"                  , emit: versions

    script:
    def VERSION = "1.1.1" // from perl script, unlikely to be updated
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $config | sed "s,hsindex = .*,hsindex = ${hisat2index}/reference," | sed "s,^faindex = .*,faindex = $genomefai," | sed "s,lastindex = .*,lastindex = ${lastindex}/reference," | sed "s,mtfaindex = .*,mtfaindex = $mtfai," | sed "s,gsize = .*,gsize = $chrsizes," | sed "s,MT_fasta = .*,MT_fasta = $mtfasta," | sed "s,refchr = .*,refchr = ${mito_name}," > new-${config}
    mkdir -p log indel bam tab bw plot
    MitoSAlt1.1.1.pl new-${config} $reads $prefix
    mv indel/*.breakpoint ${prefix}.breakpoint
    mv indel/*.cluster ${prefix}.cluster

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitosalt: $VERSION
    END_VERSIONS
    """
}
