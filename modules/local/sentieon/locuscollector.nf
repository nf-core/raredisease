process SENTIEON_LOCUSCOLLECTOR {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path('*score.txt.gz')          , emit: score        , optional: true
    tuple val(meta), path('*score.txt.gz.tbi')      , emit: score_idx    , optional: true
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def input  = bam ? '-i ' + bam.sort().join(' -i ') : ''
    """
    sentieon \\
        driver \\
        --algo LocusCollector \\
        --fun score_info ${idSample}_score.gz \\
        -t $task.cpus \\
        $input \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.score.txt.gz
    touch ${prefix}.score.txt.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
