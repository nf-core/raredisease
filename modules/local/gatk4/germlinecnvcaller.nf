process GATK4_GERMLINECNVCALLER {
    tag "$meta.id"
    label 'process_medium'
    label 'gatk4'

    input:
    tuple val(meta), path(tsv), path(ploidy)
    path  model
    val   mode

    output:
    tuple val(meta), path("*.tar"), emit: tar
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def model_command = model ? "--model $model" : ""
    def mode_command = mode ? "--run-mode $mode" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK GermlineCNVCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    tar -xvf ploidy.tar
    mkdir cnv_calls/
    gatk --java-options "-Xmx${avail_mem}g" GermlineCNVCaller \\
        $mode_command \\
        $model_command \\
        --input $tsv \\
        --contig-ploidy-calls ploidy/ \\
        --output cnv_calls/ \\
        --output-prefix $prefix \\
        $args
    tar -cvf cnv_calls.tar cnv_calls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
