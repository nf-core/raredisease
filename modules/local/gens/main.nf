process GENS {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/clinicalgenomics/gens_preproc:1.0.11'

    input:
    tuple val(meta), path(read_counts)
    tuple val(meta2), path(gvcf)
    path  gnomad_positions

    output:
    tuple val(meta), path('*.cov.bed.gz')    , emit: cov
    tuple val(meta), path('*.cov.bed.gz.tbi'), emit: cov_index
    tuple val(meta), path('*.baf.bed.gz')    , emit: baf
    tuple val(meta), path('*.baf.bed.gz.tbi'), emit: baf_index
    path  "versions.yml"                     , emit: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "The gens pre-processing module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    generate_gens_data.pl \\
        $read_counts \\
        $gvcf \\
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
