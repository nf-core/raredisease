process SENTIEON_BQSR {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

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
    def args2  = task.ext.args2 ?: ''
    def args3  = task.ext.args3 ?: ''
    def input  = bam.sort().collect{"-i $it"}.join(' ')
    def dbsnp  = known_dbsnp  ? "-k $known_dbsnp" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ ! -n \${SENTIEON_LICENSE_BASE64+x} ]; then
        echo "Initializing SENTIEON_LICENSE env variable"
        source sentieon_init.sh \${SENTIEON_LICENSE_BASE64}
    fi

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
        $args2 \\
        $input \\
        -q ${prefix}_recal_data.table \\
        --algo QualCal \\
        $dbsnp \\
        ${prefix}_recal_data.table_post \\
        --algo ReadWriter ${prefix}_recal.bam

    sentieon driver \\
        -t ${task.cpus} \\
        $args3 \\
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
