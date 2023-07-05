process SENTIEON_READWRITER {
    tag "$meta.id"
    label 'process_medium'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path('*.bam')                   , emit: bam
    tuple val(meta), path('*.bam.bai')               , emit: bai
    tuple val(meta), path('*.bam'), path('*.bam.bai'), emit: bam_bai
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def input  = bam.sort().collect{"-i $it"}.join(' ')
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
        --algo ReadWriter \\
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
