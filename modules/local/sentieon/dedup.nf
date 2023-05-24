process SENTIEON_DEDUP {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    input:
    tuple val(meta), path(bam), path(bai), path(score), path(score_idx)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path('*.bam')        , emit: bam
    tuple val(meta), path('*.bam.bai')    , emit: bai
    tuple val(meta), path('*_metrics.txt'), emit: metrics_dedup
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input  = bam.sort().collect{"-i $it"}.join(' ')
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
        --metrics ${prefix}_metrics.txt \\
        ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    touch ${prefix}_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
