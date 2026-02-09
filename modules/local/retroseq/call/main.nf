process RETROSEQ_CALL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73188352f5d5762376ee86bb739902750cc7197398f09c1b6a2b8fe3d71e22fc/data':
        'community.wave.seqera.io/library/perl-retroseq:1.5--e825fb294f7eb523' }"


    input:
    tuple val(meta), path(tab), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}"), val('retroseq'), val("1.5"), topic: versions, emit: versions_retroseq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    retroseq.pl \\
        -call \\
        $args \\
        -bam $bam \\
        -input $tab \\
        -ref $fasta \\
        -output ${prefix}.vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    """
}
