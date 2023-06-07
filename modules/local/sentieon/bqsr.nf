process SENTIEON_BQSR {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(known_dbsnp)
    tuple val(meta5), path(known_dbsnp_tbi)

    output:
    tuple val(meta), path('*.bam')       , emit: bam
    tuple val(meta), path('*.bam.bai')   , emit: bai
    tuple val(meta), path('*.table')     , emit: recal_pre
    tuple val(meta), path('*.table_post'), emit: recal_post
    tuple val(meta), path('*.csv')       , emit: recal_csv
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def args2  = task.ext.args2  ?: ''
    def args3  = task.ext.args3  ?: ''
    def dbsnp  = known_dbsnp     ? "-k $known_dbsnp" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input  = bam.sort().collect{"-i $it"}.join(' ')
    """
    if [ \${SENTIEON_LICENSE_BASE64:-"unset"} != "unset" ]; then
        echo "Initializing SENTIEON_LICENSE env variable"
        source sentieon_init.sh SENTIEON_LICENSE_BASE64
    fi

    sentieon driver  \\
        -t ${task.cpus} \\
        -r $fasta \\
        $args \\
        $input \\
        --algo QualCal \\
        $dbsnp \\
        ${prefix}.table

    sentieon driver \\
        -t ${task.cpus} \\
        -r $fasta \\
        $args2 \\
        $input \\
        -q ${prefix}.table \\
        --algo QualCal \\
        $dbsnp \\
        ${prefix}.table_post \\
        --algo ReadWriter ${prefix}.bam

    sentieon driver \\
        -t ${task.cpus} \\
        $args3 \\
        --algo QualCal \\
        --plot \\
        --before ${prefix}.table \\
        --after ${prefix}.table_post \\
        ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    touch ${prefix}.table
    touch ${prefix}.table_post
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
