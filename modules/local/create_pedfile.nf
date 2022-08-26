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

    script:
    def pedinfo = ['#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype'].join('\t')
    for(int i = 0; i<samples.size(); i++) {
        pedinfo += '\n'
        pedinfo += [samples[i].case_id, samples[i].id, samples[i].paternal, samples[i].maternal, samples[i].gender, samples[i].phenotype].join('\t');
    }
    """
    echo "$pedinfo" > family.ped
    """
}
