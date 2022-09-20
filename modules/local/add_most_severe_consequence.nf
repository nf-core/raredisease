process ADD_MOST_SEVERE_CSQ {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(vcf)
    path (variant_consequences)

    output:
    tuple val(meta), path("*.vcfparser.vcf")    , emit: vcf
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 add_most_severe_consequence.py ${vcf} ${prefix}.vcfparser.vcf ${variant_consequences}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_most_severe_consequence: v1.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcfparser.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_most_severe_consequence: v1.0
    END_VERSIONS
    """
}
