process GATK4_DETERMINEGERMLINECONTIGPLOIDY {
    tag "$meta.id"
    label 'process_medium'
    label 'gatk4'

    input:
    tuple val(meta), path(tsv)
    path  model
    // path priors

    output:
    tuple val(meta), path("ploidy.tar"), emit: tar
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def model_command = model ? "--model $model" : ""
    // def priors_command = priors ? "--contig-ploidy-priors $priors" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK DetermineGermlineContigPloidy] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    // $priors_command
    """
    gatk --java-options "-Xmx${avail_mem}g" DetermineGermlineContigPloidy \\
        $model_command \\
        --input $tsv \\
        -imr OVERLAPPING_ONLY \\
        --output ploidy/ \\
        --output-prefix $prefix \\
        $args
    tar -cvf ploidy.tar ploidy/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
