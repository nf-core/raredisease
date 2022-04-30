process SENTIEON_BQSR {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path known_dbsnp
    path known_dbsnp_tbi

    output:
    tuple val(meta), path('*_recal.bam')             , emit: bam
    tuple val(meta), path('*_recal.bam.bai')         , emit: bai
    tuple val(meta), path('*_recal_data.table')      , emit: recal_pre
    tuple val(meta), path('*_recal_data.table_post') , emit: recal_post
    tuple val(meta), path('*_recal_result.csv')      , emit: recal_csv
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def input  = bam.sort().collect{"-i $it"}.join(' ')
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
        ${prefix}_recal_data.table

    sentieon driver \\
        -t ${task.cpus} \\
        -r $fasta \\
        $args \\
        $input \\
        -q ${prefix}_recal_data.table \\
        --algo QualCal \\
        $dbsnp \\
        ${prefix}_recal_data.table_post \\
        --algo ReadWriter ${prefix}_recal.bam

    sentieon driver \\
        -t ${task.cpus} \\
        --algo QualCal \\
        --plot \\
        --before ${prefix}_recal_data.table \\
        --after ${prefix}_recal_data.table_post \\
        ${prefix}_recal_result.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_recal.bam
    touch ${prefix}_recal.bam.bai
    touch ${prefix}_recal_data.table
    touch ${prefix}_recal_data.table_post
    touch ${prefix}_recal_result.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
