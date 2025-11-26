process GENS {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/rannickscilifelab/gens_preproc:1.1.1'

    input:
    tuple val(meta), path(read_counts)
    tuple val(meta2), path(gvcf)
    path  gnomad_positions

    output:
    tuple val(meta), path('*.cov.bed')    , emit: cov
    tuple val(meta), path('*.baf.bed')    , emit: baf
    path  "versions.yml"                  , emit: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "The gens pre-processing module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    generate_gens_data.py \\
        --coverage $read_counts \\
        --gvcf $gvcf \\
        --label $prefix \\
        --baf_positions $gnomad_positions \\
        --outdir \$PWD

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_gens_data.py: \$(echo \$(generate_gens_data.py --version))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.baf.bed
    touch ${prefix}.cov.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_gens_data.py: \$(echo \$(generate_gens_data.py --version))
    END_VERSIONS
    """
}
