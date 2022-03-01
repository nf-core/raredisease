process GENS {
    tag "$meta.id"
    label 'process_medium'

    //conda (params.enable_conda ? "bioconda::gatk4=4.2.4.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'raysloks/gens_preproc:1.0.0' :
        'raysloks/gens_preproc:1.0.0' }"

    input:
    tuple val(meta), path(read_counts)
    path vcf
    path gnomad_positions

    output:
    tuple val(meta), path('*.cov.bed.gz'), emit: cov
    tuple val(meta), path('*.baf.bed.gz'), emit: baf
    path  "versions.yml"             , emit: versions

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
        $vcf \\
        $read_counts \\
        $prefix \\
        $gnomad_positions

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_gens_data.pl: \$(echo \$(generate_gens_data.pl --version 2>&1))
    END_VERSIONS
    """
}
