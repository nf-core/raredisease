process SENTIEON_DNASCOPE {
    tag	"$meta.id"
    label 'process_high'
    label 'sentieon'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path known_dbsnp
    path known_dbsnp_tbi
    path ml_model
    
    output:
    tuple val(meta), path("*_dnascope.vcf")     , emit: vcf
    tuple val(meta), path("*_dnascope.vcf.idx") , emit: vcf_index
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def dbsnp = known_dbsnp ? "-d ${known_dbsnp}" : ''
    def model = ml_model ? "--model ${ml_model}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sentieon driver \\
        -t $task.cpus \\
        -r $fasta \\
        $args \\
        -i $bam \\
        --algo DNAscope \\
        $dbsnp \\
        $args2 \\
        $args3 \\
        $model \\
        ${prefix}_dnascope.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_dnascope.vcf
    touch ${prefix}_dnascope.vcf.idx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g" )
    END_VERSIONS
    """
}
