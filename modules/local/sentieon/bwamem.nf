process SENTIEON_BWAMEM {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    input:
    tuple val(meta), path(reads)
    path fasta
    path fai
    tuple val(meta2), path(index) // meta2 has same purpose as meta, and holds information about the genome/index

    output:
    tuple val(meta), path('*.bam'), emit: bam
    tuple val(meta), path('*.bai'), emit: bai
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    if [ \${SENTIEON_LICENSE_BASE64:-"unset"} != "unset" ]; then
        echo "Initializing SENTIEON_LICENSE env variable"
        source sentieon_init.sh SENTIEON_LICENSE_BASE64
    fi

    sentieon bwa mem \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        $args \\
        | sentieon \\
            util \\
            sort \\
            -r $fasta \\
            -o ${prefix}.bam \\
            -t $task.cpus \\
            $args2 \\
            --sam2bam \\
            -i -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$(sentieon bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
