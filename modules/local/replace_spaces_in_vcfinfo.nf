process REPLACE_SPACES_IN_VCFINFO {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*_reformatted.vcf"), emit: vcf
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 <<CODE
    from pathlib import Path
    with open("${input}") as input:
        output_fn = Path("${input}").stem + "_reformatted.vcf"
        with open(output_fn,'w') as output:
            for line in input:
                if line.startswith("#"):
                    output.write(line)
                else:
                    spl = line.strip().split("\\t")
                    spl[7] = spl[7].replace(" ","_")
                    output.write(("\\t").join(spl)+"\\n")
    CODE

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        replace_spaces_in_vcfinfo: v1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    python3 <<CODE
    from pathlib import Path
    with open("${input}") as input:
        output_fn = Path("${input}").stem + "_reformatted.vcf"
        with open(output_fn,'w') as output:
            pass
    CODE

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        replace_spaces_in_vcfinfo: v1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
