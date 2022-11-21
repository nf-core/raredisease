process CHECK_INPUT_VCF {
    tag "check_vcf"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    path vcf

    output:
    path '*.csv'         , emit: csv
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/raredisease/bin/
    """
    export INPUT_FILE=${vcf}
    export OUTPUT_FILE="checked_vcfs.csv"

    python3 <<CODE
    import os, gzip
    file_in  = os.environ.get('INPUT_FILE')
    file_out = os.environ.get('OUTPUT_FILE')
    if file_in.endswith(".gz"):
        with open(file_out,'w') as out:
            base = os.path.basename(file_in).rsplit(".",2)[0]
            out.write("id,filepath,processed\\n")
            with gzip.open(file_in,'rt') as vcf:
                for line in vcf:
                    if line.startswith("##bcftools_norm"):
                        out.write(base + "," + os.path.abspath(file_in) + ",yes\\n")
                        break
                    elif not line.startswith("#"):
                        out.write(base + "," + os.path.abspath(file_in) + ",no\\n")
                        break
    else:
        print("Please compress %s using bgzip" %file_in)
    CODE

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        change_input_vcf: v1.0
    END_VERSIONS
    """

    stub:
    """
    touch checked_vcfs.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        check_input_vcf: v1.0
    END_VERSIONS
    """
}
