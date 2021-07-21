// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CALCULATE_MAPPED {
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"samtools", meta:[:], publish_by_meta:['id']) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    tuple val(meta), path(stats)

    output:
    path "*.stats"     , emit: stats

    script:
    def prefix = stats.getSimpleName()
    """
    cat $stats | grep -E "raw total sequences:|reads mapped:" | cut -f3 | tr '\\n' '\\t' | awk '{print (\$2/\$1) * 100}' | while read line; do awk -v p=\$line 'NR==15 {print "percentage mapped reads:\\t" p }1' $stats; done > ${prefix}.stats
    """
}
