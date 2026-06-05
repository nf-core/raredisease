process FILTERVEP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/48/481e0a88c723eb6abe5a6c36e50bc2119b460f3c6f71aff2111a22e8b704c5d9/data' :
        'community.wave.seqera.io/library/cyvcf2:0.32.1--569b36b775b7f1e5' }"

    input:
    tuple val(meta), path(input)
    path feature_file

    output:
    tuple val(meta), path("*.${extension}"), emit: output
    tuple val("${task.process}"), val('filtervep'), val('1.0.0'), topic: versions, emit: versions_filtervep

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    extension     = args.contains("--format tab") ? "txt" : "vcf"
    """
    filtervep.py \\
        ${args} \\
        --input_file ${input} \\
        --output_file ${prefix}.${extension} \\
        --only_matched
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    extension  = args.contains("--format tab") ? "txt" : "vcf"
    """
    touch ${prefix}.${extension}
    """
}
