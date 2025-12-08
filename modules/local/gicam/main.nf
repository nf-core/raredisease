process GICAM_INFER {
    // https://github.com/Clinical-Genomics/rdds/tree/master/src/rdds/gicam

    tag "${meta.id}"
    label 'process_high'

    container "docker.io/clinicalgenomics/rdds_mivmir:v1.12.0-rc6"

    beforeScript "mkdir ${task.workDir}/rdds-tmp"
    afterScript "rm -r ${task.workDir}/rdds-tmp"
    containerOptions {[
        workflow.containerEngine.equals("singularity") ? "--bind ${task.workDir}/rdds-tmp:/rdds/tmp" : "",
        workflow.containerEngine.equals("docker") ? "--tmpfs /rdds/tmp": "",
        ""
    ].minus("").join(" ")}

    input:
    tuple val(meta), path(input_vcf)

    output:
    tuple val(meta), path('*-predictions.vcf'), emit: vcf
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = 'v1.12.0-rc6'
    """
    . /opt/pyenv/bin/activate
    export PYTHONPATH=/rdds/src
    python3 -m rdds.gicam infer-vcf --cpu_cores ${task.cpus} ${input_vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gicam: ${VERSION}
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
