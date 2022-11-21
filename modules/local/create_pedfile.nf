process MAKE_PED {
    tag "make_ped"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    val(samples)

    output:
    path '*.ped'          , emit: ped
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    exec:
    def pedfile  = task.workDir.resolve('family.ped')
    pedfile.text = ['#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype'].join('\t')
    def samples_list = []
    for(int i = 0; i<samples.size(); i++) {
        sample_tokenized   =  samples[i].id.tokenize("_")
        sample_tokenized.removeLast()
        sample_name        =  sample_tokenized.join("_")
        if (!samples_list.contains(sample_name)) {
            pedfile.append('\n' + [samples[i].case_id, sample_name, samples[i].paternal, samples[i].maternal, samples[i].gender, samples[i].phenotype].join('\t'));
            samples_list.add(sample_name)
        }
    }
    def versionfile = task.workDir.resolve('versions.yml')
    versionfile.text = '"' + task.process + '"' + '\n    create_pedfile: 1.0\n'
}
