process ENSEMBLVEP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::ensembl-vep=105.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:105.0--pl5321h4a94de4_1' :
        'quay.io/biocontainers/ensembl-vep:105.0--pl5321h4a94de4_1' }"

    input:
    tuple val(meta), path(vcf), path(index)
    val   genome
    val   species
    val   cache_version
    path  cache

    output:
    tuple val(meta), path("*.ann.vcf.gz"), emit: vcf
    path "*.summary.html"                , emit: report, optional: true
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"
    """
    mkdir $prefix

    vep \\
        -i $vcf \\
        -o ${prefix}.ann.vcf.gz \\
        $args \\
        --assembly $genome \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache $dir_cache \\
        --fork $task.cpus \\
        --vcf \\
        --compress_output bgzip \\
        --stats_file ${prefix}.summary.html

    rm -rf $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ann.vcf.gz
    touch ${prefix}.summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
