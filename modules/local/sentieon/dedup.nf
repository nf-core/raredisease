process SENTIEON_DEDUP {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    input:
    tuple val(meta), path(bam), path(bai), path(score), path(score_idx)
    path fasta
    path fai

    output:
    tuple val(meta), path('*dedup.bam')          , emit: bam
    tuple val(meta), path('*dedup.bam.bai')      , emit: bai
    tuple val(meta), path('*dedup_metrics.txt')  , emit: metrics_dedup
    path  "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def input = bam.sort().collect{"-i $it"}.join(' ')
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ \${SENTIEON_LICENSE_BASE64:-"unset"} != "unset" ]; then
        echo "Initializing SENTIEON_LICENSE env variable"
        source sentieon_init.sh SENTIEON_LICENSE_BASE64
    fi

    sentieon \\
        driver \\
        -t $task.cpus \\
        $input \\
        $args \\
        --algo Dedup \\
        --score_info $score \\
        --metrics ${prefix}_dedup_metrics.txt \\
        ${prefix}_dedup.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_dedup.bam
    touch ${prefix}_dedup.bam.bai
    touch ${prefix}_dedup_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
