process SENTIEON_BQSR {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    input:
    tuple val(meta), path(bam), path(bai), path(score), path(score_idx)
    path fasta
    path fai
    path known_dbsnp

    output:
    tuple val(meta), path('*.recal.bam')            , emit: bam
    tuple val(meta), path('*.recal.bai')            , emit: bai
    tuple val(meta), path('*recal_data.table')      , emit: recal_pre
    tuple val(meta), path('*recal_data.table.post') , emit: recal_post
    tuple val(meta), path('*recal_result.csv')      , emit: recal_csv
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def input  = bam ? '-i ' + bam.sort().join(' -i ') : ''
    def dbsnp  = known_dbsnp  ? "-k $known_dbsnp" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sentieon driver  \\
        -t ${task.cpus} \\
        -r $fasta \\
        $args \\
        $input \\
        --algo QualCal \\
        $dbsnp \\
        ${prefix}.recal_table.table

    sentieon driver \\
        -t ${task.cpus} \\
        -r $fasta \\
        $args \\
        $input \\
        -q ${prefix}.recal_data.table \\
        --algo QualCal \\
        $dbsnp \\
        $prefix.recal_data.table.post \\
        --algo ReadWriter ${prefix}.recal.bam

    sentieon driver \\
        -t ${task.cpus} \\
        --algo QualCal \\
        --plot \\
        --before ${prefix}.recal.table \\
        --after ${prefix}.table.post \\
        ${prefix}_recal_result.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.recal.bam
    touch ${prefix}.recal.bai
    touch ${prefix}.recal_data.table
    touch ${prefix}.recal_data.table.post
    touch ${prefix}_recal_result.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
