process ENSEMBLVEP_VEP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::ensembl-vep=110.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f79e4fcf077482a6f17596bd68f7ba5a8c1972a439a5b7899549bc33baed8c9d/data' :
        'community.wave.seqera.io/library/ensembl-vep:110.0--cf66340b49c82813' }"

    input:
    tuple val(meta), path(vcf), path(custom_extra_files)
    val   genome
    val   species
    val   cache_version
    path  cache
    tuple val(meta2), path(fasta)
    path  extra_files

    output:
    tuple val(meta), path("*.vcf.gz")  , optional:true, emit: vcf
    tuple val(meta), path("*.tab.gz")  , optional:true, emit: tab
    tuple val(meta), path("*.json.gz") , optional:true, emit: json
    tuple val(meta), val("${task.process}"), val('ensemblvep'), path("*.html"), topic: multiqc_files, emit: report, optional: true
    tuple val("${task.process}"), val('ensemblvep'), eval("vep --help | sed -n '/ensembl-vep/s/.*: //p'"), topic: versions, emit: versions_ensemblvep

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json")? 'json' : args.contains("--tab")? 'tab' : 'vcf'
    def compress_cmd = args.contains("--compress_output") ? '' : '--compress_output bgzip'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"
    def reference = fasta ? "--fasta $fasta" : ""
    """
    vep \\
        -i $vcf \\
        -o ${prefix}.${file_extension}.gz \\
        $args \\
        $compress_cmd \\
        $reference \\
        --assembly $genome \\
        --species $species \\
        --cache \\
        --cache_version $cache_version \\
        --dir_cache $dir_cache \\
        --fork $task.cpus \\
        --stats_file ${prefix}.summary.html \\

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.tab.gz
    touch ${prefix}.json.gz
    touch ${prefix}.summary.html

    """
}
