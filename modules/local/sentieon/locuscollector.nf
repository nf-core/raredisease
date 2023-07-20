process SENTIEON_LOCUSCOLLECTOR {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path('*txt.gz')    , emit: score    , optional: true
    tuple val(meta), path('*txt.gz.tbi'), emit: score_idx, optional: true
    path  "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def input  = bam.sort().collect{"-i $it"}.join(' ')
    def prefix = task.ext.prefix ? "${task.ext.prefix}.txt.gz" : "${meta.id}.txt.gz"
    """
    if [ \${SENTIEON_LICENSE_BASE64:-"unset"} != "unset" ]; then
        echo "Initializing SENTIEON_LICENSE env variable"
        source sentieon_init.sh SENTIEON_LICENSE_BASE64
    fi

    sentieon \\
        driver \\
        -t $task.cpus \\
        $input \\
        --algo LocusCollector \\
        --fun score_info $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ? "${task.ext.prefix}.txt.gz" : "${meta.id}.txt.gz"
    """
    touch ${prefix}.txt.gz
    touch ${prefix}.txt.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
