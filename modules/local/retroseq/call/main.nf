process RETROSEQ_CALL {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::perl-retroseq=1.5=pl5321hdfd78af_1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker.io/clinicalgenomics/retroseq:1.5_9d4f3b5-1' : 'docker.io/clinicalgenomics/retroseq:1.5_9d4f3b5-1' }"


    input:
    tuple val(meta), path(tab), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.5"

    """
    retroseq.pl \\
        -call \\
        $args \\
        -bam $bam \\
        -input $tab \\
        -ref $fasta \\
        -output ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        retroseq_call: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.5"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        retroseq_call: $VERSION
    END_VERSIONS
    """
}
