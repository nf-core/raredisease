process SPLIT_CHR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--h5bf99c6_6':
        'biocontainers/bioawk:1.0--h5bf99c6_6' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("chr_split"), emit: output
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '1.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir -p chr_split
    bioawk -v outdir="chr_split" 'BEGIN{RS=">"; FS="\\n"} NR>1 {fnme=outdir "/" \$1 ".fa"; print ">" \$0 > fnme; close(fnme); if (NR==25) exit;}' ${input}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioawk: $VERSION
    END_VERSIONS
    """
}
