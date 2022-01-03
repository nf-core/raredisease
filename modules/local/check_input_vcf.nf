process CHECK_INPUT_VCF {
    tag "check_vcf"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

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
