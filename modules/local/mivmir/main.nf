process MIVMIR_INFER {
    // https://github.com/Clinical-Genomics/rdds/tree/master/src/rdds/variant_rank_score

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
    tuple val("${task.process}"), val('mivmir'), val('v1.12.0-rc6'), topic: versions, emit:versions_mivmir

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    . /opt/pyenv/bin/activate
    export PYTHONPATH=/rdds/src
    python3 -m rdds.variant_rank_score predict-on-vcf  --cpu_cores ${task.cpus} ${input_vcf}
    """
}

process MIVMIR_INTERNAL_UNIT_TEST {
    // Test inference API and numerical reproducibility
    container "docker.io/clinicalgenomics/rdds_mivmir:v1.12.0-rc6"

    beforeScript "mkdir ${task.workDir}/rdds-tmp"
    afterScript "rm -r ${task.workDir}/rdds-tmp"
    containerOptions {[
        workflow.containerEngine.equals("singularity") ? "--bind ${task.workDir}/rdds-tmp:/rdds/tmp" : "",
        workflow.containerEngine.equals("docker") ? "--tmpfs /rdds/tmp": "",
        ""
    ].minus("").join(" ")}

    script:
    """
    . /opt/pyenv/bin/activate
    export PYTHONPATH=/rdds/src
    export TF_ENABLE_ONEDNN_OPTS=0
    cd /rdds/src/tests
    python3 -m pytest -x -v -s -o log_cli=true --full-trace lib/determinism variant_rank_score/test_inference.py
    """
}
