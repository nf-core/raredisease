// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GLNEXUS_MERGE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

//    conda (params.enable_conda ? "YOUR-TOOL-HERE" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://clinicalgenomics/glnexus:v1.4.1"
    } else {
        container "clinicalgenomics/glnexus:v1.4.1"
    }

    input:
    tuple val(meta), path(gvcfs)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    glnexus_cli \\
        --threads $task.cpus --mem-gbytes $task.memory \\
        $options.args \\
        $gvcfs \\
        > ${prefix}.bcf 

    echo \$(glnexus_cli) | | head -n 1 | sed 's/^.*release //; s/ .*$//' > ${software}.version.txt
    """
}
