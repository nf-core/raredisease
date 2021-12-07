// Import generic module functions
// include { saveFiles } from './functions'

// params.options = [:]

process CHECK_INPUT_VCF {
    tag "check_vcf"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'genome', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path vcf

    output:
    path '*.txt'       , emit: txt

    script: // This script is bundled with the pipeline, in nf-core/raredisease/bin/
    """
    check_input_vcf.py \\
        --INPUT_VCF $vcf \\
        --OUTPUT checked_vcfs.txt
    """
}
