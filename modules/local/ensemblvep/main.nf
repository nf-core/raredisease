process ENSEMBLVEP {
    tag "$meta.id"
    label 'process_medium'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local VEP module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    container "docker.io/ensemblorg/ensembl-vep:release_107.0"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta)
    val   genome
    val   species
    val   cache_version
    path  cache
    path  extra_files

    output:
    tuple val(meta), path("*.vcf")    , optional:true, emit: vcf
    tuple val(meta), path("*.tab")    , optional:true, emit: tab
    tuple val(meta), path("*.json")   , optional:true, emit: json
    tuple val(meta), path("*.vcf.gz") , optional:true, emit: vcf_gz
    tuple val(meta), path("*.tab.gz") , optional:true, emit: tab_gz
    tuple val(meta), path("*.json.gz"), optional:true, emit: json_gz
    path "*.summary.html"             , optional:true, emit: report
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json")? 'json' : args.contains("--tab")? 'tab' : 'vcf'
    def compress_out = args.contains("--compress_output") ? '.gz' : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stats_file     = args.contains("--no_stats") ? '' : "--stats_file ${prefix}.summary.html"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"
    def reference = fasta ? "--fasta $fasta" : ""

    """
    vep \\
        -i $vcf \\
        -o ${prefix}.${file_extension}${compress_out} \\
        $args \\
        $reference \\
        --assembly $genome \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache $dir_cache \\
        --fork $task.cpus \\
        ${stats_file}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    touch ${prefix}.tab
    touch ${prefix}.json
    touch ${prefix}.vcf.gz
    touch ${prefix}.tab.gz
    touch ${prefix}.json.gz
    touch ${prefix}.summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
