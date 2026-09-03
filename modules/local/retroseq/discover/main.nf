process RETROSEQ_DISCOVER {
    tag "$meta.id"
    label 'process_low'

    container 'docker.io/clinicalgenomics/retroseq:1.5_9d4f3b5-1'

    input:
    tuple val(meta), path(bam), path(bai)
    path(me_references)
    val(me_types)

    output:
    tuple val(meta), path("*.tab"), emit: tab
    tuple val("${task.process}"), val('retroseq'), val("1.5"), topic: versions, emit: versions_retroseq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    paste <(printf "%s\\n" $me_types | tr -d '[],') <(printf "%s\\n" $me_references) > me_reference_manifest.tsv
    retroseq.pl \\
        -discover \\
        $args \\
        -bam $bam \\
        -refTEs me_reference_manifest.tsv\\
        -output ${prefix}.tab

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    paste <(printf "%s\\n" $me_types | tr -d '[],') <(printf "%s\\n" $me_references) > me_reference_manifest.tsv
    touch ${prefix}.tab
    """
}
