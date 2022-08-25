process MAKE_PED {
    tag "make_ped"

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    val(samples)

    output:
    path '*.ped'       , emit: ped

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/raredisease/bin/
    def pedinfo = ""
    for(int i = 0; i<samples.size(); i++) {
        pedinfo += samples[i].case_id + '\t' + samples[i].id + '\t' + samples[i].paternal + '\t' + samples[i].maternal + '\t' + samples[i].gender + '\t' + samples[i].phenotype + '\n';
    }
    """
    echo "#family_id\tsample_id\tfather\tmother\tsex\tphenotype" >family.ped
    echo "$pedinfo" >>family.ped
    """
}
