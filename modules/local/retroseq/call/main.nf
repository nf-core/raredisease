process RETROSEQ_CALL {
    tag "$meta.id"
    label 'process_low'

    container 'docker.io/clinicalgenomics/retroseq:1.5_9d4f3b5-1'

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
