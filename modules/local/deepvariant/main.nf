process DEEPVARIANT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::deepvariant=1.2.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://google/deepvariant:1.2.0' :
        'google/deepvariant:1.2.0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    tuple val(fai_meta), path(fasta_fai)

    output:
    tuple val(meta), path("*.vcf.gz"),  emit: vcf
    tuple val(meta), path("*g.vcf.gz"), emit: gvcf
    path  "versions.yml",               emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    /opt/deepvariant/bin/run_deepvariant \\
        --ref=$fasta \\
        --reads=$bam \\
        --output_vcf=${prefix}.vcf.gz \\
        --output_gvcf=${prefix}.g.vcf.gz \\
        $args \\
        --num_shards=${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """

}
