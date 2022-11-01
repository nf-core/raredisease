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
    path (pli_gene)

    output:
    tuple val(meta), path("*.csq_pli.vcf")    , emit: csq_pli_vcf
    tuple val(meta), path("*.csq.vcf")        , emit: csq_vcf
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 add_most_severe_consequence.py --file_in ${vcf} --file_out ${prefix}.csq.vcf --variant_csq ${variant_consequences}
    python3 add_most_severe_pli.py --file_in ${prefix}.csq.vcf --file_out ${prefix}.csq_pli.vcf --pli ${pli}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_most_severe_consequence_pli: v1.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csq_pli.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        add_most_severe_consequence_pli: v1.0
    END_VERSIONS
    """
}
