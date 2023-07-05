process SENTIEON_DNASCOPE {
    tag	"$meta.id"
    label 'process_high'
    label 'sentieon'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(known_dbsnp)
    tuple val(meta5), path(known_dbsnp_tbi)
    path call_interval
    path ml_model

    output:
    tuple val(meta), path("*.vcf.gz")                      , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi")                  , emit: index
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf_index
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args   ?: ''
    def args2    = task.ext.args2  ?: ''
    def interval = call_interval   ? "--interval ${call_interval}" : ''
    def dbsnp    = known_dbsnp     ? "-d ${known_dbsnp}"           : ''
    def model    = ml_model        ? "--model ${ml_model}"         : ''
    def prefix   = task.ext.prefix ?: "${meta.id}"

    """
    sentieon driver \\
        -t $task.cpus \\
        -r $fasta \\
        $args \\
        -i $bam \\
        --algo DNAscope \\
        $dbsnp \\
        $interval \\
        $args2 \\
        $model \\
        ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g" )
    END_VERSIONS
    """
}
