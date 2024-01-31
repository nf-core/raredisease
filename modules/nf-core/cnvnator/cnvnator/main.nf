process CNVNATOR_CNVNATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvnator:0.4.1--py310h2dce045_7':
        'biocontainers/cnvnator:0.4.1--py310h2dce045_7' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(root)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)
    val step

    output:
    tuple val(output_meta), path("${output_meta.id}_${step}.root"), emit: root
    tuple val(output_meta), path("*.tab")                         , emit: tab, optional: true
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def input_cmd = bam             ? "-tree ${bam}"      : ''
    output_meta   = bam             ? meta                : meta2
    prefix        = task.ext.prefix ?: bam ? "${meta.id}" : "${meta2.id}"
    if (fasta) {
        reference = fasta.isDirectory() ? "-d ${fasta}" : "-fasta ${fasta}"
    } else {
        reference = ''
    }
    calls_cmd = args.contains("-call") ? "> ${prefix}.tab" : ''
    mv_cmd    = "mv ${prefix}.root ${prefix}_${step}.root"
    steps     = ["his", "stat", "partition", "call"]
    cp_cmd    = steps.contains(step) ? "cp ${root} ${prefix}.root" :""
    """
    $cp_cmd
    cnvnator \\
        -root ${prefix}.root \\
        $args \\
        $reference \\
        $input_cmd \\
        $calls_cmd

    $mv_cmd 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CNVnator: \$(echo \$(cnvnator 2>&1 | sed -n '3p' | sed 's/CNVnator v//'))
    END_VERSIONS
    """

    stub:
    def args      = task.ext.args   ?: ''
    prefix        = task.ext.prefix ?: bam ? "${meta.id}" : "${meta2.id}"
    output_meta   = bam             ? meta                : meta2
    def calls_cmd = args.contains("-call") ? "touch ${prefix}.tab" : ''
    """
    touch ${prefix}_${step}.root
    $calls_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CNVnator: \$(echo \$(cnvnator 2>&1 | sed -n '3p' | sed 's/CNVnator v//'))
    END_VERSIONS
    """
}
