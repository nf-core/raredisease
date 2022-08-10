process FILTER_VEP {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::ensembl-vep=105.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:105.0--pl5321h4a94de4_1' :
        'quay.io/biocontainers/ensembl-vep:105.0--pl5321h4a94de4_1' }"

    input:
    tuple val(meta), path(vcf)
    path (select_feature_file_bed)

    output:
    tuple val(meta), path("*.ann_filter.vcf.gz"), emit: vcf
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export INPUT_FILE=${select_feature_file_bed}
    export OUTPUT_FILE="feature.list"

    python3 <<CODE
    import os
    file_in  = os.environ.get('INPUT_FILE')
    file_out = os.environ.get('OUTPUT_FILE')
    with open(file_out,'w') as out:
        with open(file_in,'rt') as bed:
            for line in bed:
                if line.startswith("#"):
                    pass
                else:
                    columns = line.strip().split(",")
                    out.write(columns[3] + "\\n")
    CODE

    filter_vep \\
        --format vcf \\
        --input_file $vcf \\
        --output_file ${prefix}.ann_filter.vcf.gz \\
        --filter \"SYMBOL in feature.list\" \\
        --soft_filter

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ann_filter.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
