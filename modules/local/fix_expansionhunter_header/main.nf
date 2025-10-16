process BCFTOOLS_FIXEXPANSIONHUNTERHEADER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create header file with missing INFO fields from ExpansionHunter
    cat > expansionhunter_headers.txt << 'EOF'
    ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
    ##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference copy number">
    ##INFO=<ID=REPID,Number=1,Type=String,Description="Repeat identifier from the variant catalog">
    ##INFO=<ID=RL,Number=1,Type=Integer,Description="Reference length in bp">
    ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit in the reference orientation">
    EOF

    # Add missing headers to VCF
    bcftools annotate \\
        -h expansionhunter_headers.txt \\
        $args \\
        -o ${prefix}.vcf \\
        $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
