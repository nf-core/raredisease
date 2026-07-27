process MT_DELETION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b9/b946e2f0e77ec69853787dfc8b312bd7e9d5c65a11a613ce918469a9566992e3/data' :
        'community.wave.seqera.io/library/samtools:1.23--12d9384dd0649f36' }"

    input:
    tuple val(meta), path(input), path(input_index)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path('*.txt'), emit: mt_del_result
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools
    tuple val("${task.process}"), val('mitodel'), val("1.0"), topic: versions, emit: versions_mitodel

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    samtools stats --threads ${task.cpus} $args ${reference} ${input} | \\
    grep -E ^IS | \\
    awk 'BEGIN {sum=0} (\$2>=1200 && \$2<=15000) {sum=sum+\$3} (\$2<1200 || \$2>15000) {sum_norm=sum_norm+\$3} END \\
    {print "intermediate discordant ", sum, "normal ", sum_norm, "ratio ppk", sum*1000/(sum_norm+sum)}' 1> ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
