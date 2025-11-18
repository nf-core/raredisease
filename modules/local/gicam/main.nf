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
    tuple val("${task.process}"), val('gicam'), val('v1.12.0-rc6'), topic: versions, emit: versions_gicam

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    . /opt/pyenv/bin/activate
    export PYTHONPATH=/rdds/src
    python3 -m rdds.gicam infer-vcf --cpu_cores ${task.cpus} ${input_vcf}
    """
}
