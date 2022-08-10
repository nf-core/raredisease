process MAKE_PED {
    tag "make_ped"

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    path samplesheet

    output:
    path '*.ped'       , emit: ped

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/raredisease/bin/
    """
    export INPUT_FILE=${samplesheet}
    export OUTPUT_FILE="familyinfo.ped"

    python3 <<CODE
    import os
    file_in  = os.environ.get('INPUT_FILE')
    file_out = os.environ.get('OUTPUT_FILE')
    sample_dict = {}
    with open(file_out,'w') as out:
        out.write("#family_id\\tsample_id\\tfather\\tmother\\tsex\\tphenotype\\n")
        infile = open(file_in).readlines()[1:]
        for line in infile:
            columns = line.strip().split(",")
            sample_dict[columns[0]] = columns[8] + "\\t" + columns[0] + "\\t" + columns[6] + "\\t" + columns[7] + "\\t" + columns[4] + "\\t" + columns[5]
        for i in sample_dict:
            out.write(sample_dict[i] + "\\n")
    CODE
    """
}
