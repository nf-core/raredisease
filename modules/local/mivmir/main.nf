process MIVMIR_INFER {
    // https://github.com/Clinical-Genomics/rdds/tree/master/src/rdds/variant_rank_score

    tag "${meta.id}"
    label 'process_high'

    container "docker.io/clinicalgenomics/rdds_mivmir:v1.14.0"

    beforeScript "mkdir ${task.workDir}/rdds-tmp"
    afterScript "rm -r ${task.workDir}/rdds-tmp"
    containerOptions {[
        workflow.containerEngine.equals("singularity") ? "--bind ${task.workDir}/rdds-tmp:/rdds/tmp" : "",
        workflow.containerEngine.equals("docker") ? "--tmpfs /rdds/tmp": "",
        ""
    ].minus("").join(" ")}

    input:
    tuple val(meta), path(input_vcf)
    val run_internal_test

    output:
    tuple val(meta), path('*-predictions.vcf'), emit: vcf
    tuple val("${task.process}"), val('mivmir'), val('v1.14.0'), topic: versions, emit:versions_mivmir

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    set -e
    . /opt/pyenv/bin/activate
    export PYTHONPATH=/rdds/src
    if [ "$run_internal_test" == "true" ]; then
        # Test inference API and numerical reproducibility
        python3 -m pytest /rdds/src/tests/variant_rank_score -k test_inference
    fi
    python3 -m rdds.variant_rank_score predict-on-vcf \
    --cpu_cores ${task.cpus} \
    $args \
    ${input_vcf}
    """
}
