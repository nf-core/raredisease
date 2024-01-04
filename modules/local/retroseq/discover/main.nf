process RETROSEQ_DISCOVER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::perl-retroseq=1.5=pl5321hdfd78af_1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker.io/clinicalgenomics/retroseq:1.5_9d4f3b5-1' : 'docker.io/clinicalgenomics/retroseq:1.5_9d4f3b5-1' }"


    input:
    tuple val(meta), path(bam), path(bai)
    path(me_references)
    val(me_types)

    output:
    tuple val(meta), path("*.tab"), emit: tab
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.5" // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.

    """
    paste <(printf "%s\\n" $me_types | tr -d '[],') <(printf "%s\\n" $me_references) > me_reference_manifest.tsv
    retroseq.pl \\
        -discover \\
        $args \\
        -bam $bam \\
        -refTEs me_reference_manifest.tsv\\
        -output ${prefix}.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        retroseq_discover: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.5"
    """
    paste <(printf "%s\\n" $me_types | tr -d '[],') <(printf "%s\\n" $me_references) > me_reference_manifest.tsv
    touch ${prefix}.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        retroseq_discover: $VERSION
    END_VERSIONS
    """
}
