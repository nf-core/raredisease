process SENTIEON_DNAMODELAPPLY {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    input:
    tuple val(meta), path(vcf), path(vcf_idx)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    path ml_model

    output:
    tuple val(meta), path("*.vcf.gz")                      , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi")                  , emit: index
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf_index
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sentieon driver \\
        -t $task.cpus \\
        -r $fasta \\
        --algo DNAModelApply \\
        --model $ml_model \\
        -v $vcf \\
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
