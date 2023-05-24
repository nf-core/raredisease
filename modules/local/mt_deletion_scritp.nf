process MT_DELETION {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path('*.txt'), emit: mt_del_script
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    grep -E ^IS $stats | awk 'BEGIN {sum=0} (\$2>=1200 && \$2<=15000) {sum=sum+\$3} (\$2<1200 || \$2>15000) {sum_norm=sum_norm+\$3} END \\
    {print "intermediate discordant ", sum, "normal ", sum_norm, "ratio ppk", sum*1000/(sum_norm+sum)}' 1> ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mt_del: v1.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mt_del.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mt_del: v1.0
    END_VERSIONS
    """
}
