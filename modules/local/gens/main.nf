process GENS {
    tag "$meta.id"
    label 'process_medium'

    container 'raysloks/gens_preproc:1.0.1'

    input:
    tuple val(meta), path(read_counts)
    path  vcf
    path  gnomad_positions

    output:
    tuple val(meta), path('*.cov.bed.gz'), emit: cov
    tuple val(meta), path('*.baf.bed.gz'), emit: baf
    path  "versions.yml"                 , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 6
    if (!task.memory) {
        log.info '[Gens] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    generate_gens_data.pl \\
        $read_counts \\
        $vcf \\
        $prefix \\
        $gnomad_positions

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_gens_data.pl: \$(echo \$(generate_gens_data.pl --version 2>&1))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.baf.bed.gz
    touch ${prefix}.baf.bed.gz.tbi
    touch ${prefix}.cov.bed.gz
    touch ${prefix}.cov.bed.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_gens_data.pl: \$(echo \$(generate_gens_data.pl --version 2>&1))
    END_VERSIONS
    """
}
