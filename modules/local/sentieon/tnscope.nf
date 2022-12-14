process SENTIEON_TNSCOPE {
    tag	"$meta.id"
    label 'process_high'
    label 'sentieon'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("*_TNscope_MTcalls_unfiltered.vcf.gz")        , emit: vcf
    tuple val(meta), path("*_TNscope_MTcalls_unfiltered.vcf.gz.tbi")    , emit: vcf_index
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def interval      = task.ext.args ?: ''
    def call_settings = task.ext.args2 ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"

    """
    sentieon driver \\
        -t $task.cpus \\
        -r $fasta \\
        -i $bam \\
        $interval \\
        --algo TNscope \\
        --tumor_sample ${meta.id} \\
        $call_settings \\
        ${prefix}_TNscope_MTcalls_unfiltered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_TNscope_MTcalls_unfiltered.vcf.gz
    touch ${prefix}_TNscope_MTcalls_unfiltered.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g" )
    END_VERSIONS
    """
}
