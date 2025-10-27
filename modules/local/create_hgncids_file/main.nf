process CREATE_HGNCIDS_FILE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(input)

    output:
    path("*_reformatted.txt"), emit: txt
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 <<CODE
    from pathlib import Path
    with open("${input}") as input:
        output_fn = Path("${input}").stem + "_reformatted.txt"
        with open(output_fn,'w') as output:
            if "scout" == "${meta.id}":
                for line in input:
                    if not line.startswith("#") and line.strip():
                        spl = line.strip().split("\\t")
                        output.write(spl[3]+"\\n")
            else:
                for line in input:
                    if not line.startswith("#"):
                        output.write(line)
    CODE

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_hgncids_file: v1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    python3 <<CODE
    from pathlib import Path
    with open("${input}") as input:
        output_fn = Path("${input}").stem + "_reformatted.txt"
        with open(output_fn,'w') as output:
            pass
    CODE

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_hgncids_file: v1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
