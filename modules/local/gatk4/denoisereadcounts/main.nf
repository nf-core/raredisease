process GATK4_DENOISEREADCOUNTS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(read_counts)
    path panel_of_normals

    output:
    tuple val(meta), path('*.standardizedCR.tsv'), emit: standardized_read_counts
    tuple val(meta), path('*.denoisedCR.tsv')    , emit: denoised_read_counts
    path  "versions.yml"                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 12288
    if (!task.memory) {
        log.info '[GATK DenoiseReadCounts] Available memory not known - defaulting to 12GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M" DenoiseReadCounts \\
        -I $read_counts \\
        --count-panel-of-normals $panel_of_normals \\
        --standardized-copy-ratios ${prefix}.standardizedCR.tsv \\
        --denoised-copy-ratios ${prefix}.denoisedCR.tsv \\
        $args \\
        --tmp-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.standardizedCR.tsv
    touch ${prefix}.denoisedCR.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
