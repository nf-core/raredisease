process RETROSEQ_DISCOVER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73188352f5d5762376ee86bb739902750cc7197398f09c1b6a2b8fe3d71e22fc/data':
        'community.wave.seqera.io/library/perl-retroseq:1.5--e825fb294f7eb523' }"


    input:
    tuple val(meta), path(bam), path(bai)
    path(me_references)
    val(me_types)

    output:
    tuple val(meta), path("*.tab"), emit: tab
    tuple val("${task.process}"), val('retroseq'), val("1.5"), topic: versions, emit: versions_retroseq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    paste <(printf "%s\\n" $me_types | tr -d '[],') <(printf "%s\\n" $me_references) > me_reference_manifest.tsv
    retroseq.pl \\
        -discover \\
        $args \\
        -bam $bam \\
        -refTEs me_reference_manifest.tsv\\
        -output ${prefix}.tab

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    paste <(printf "%s\\n" $me_types | tr -d '[],') <(printf "%s\\n" $me_references) > me_reference_manifest.tsv
    touch ${prefix}.tab
    """
}
