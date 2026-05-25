process SALTSHAKER_TO_HTML {
    tag "$meta.caseid"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a2/a23c958d5a0439419f82069df6217c542f6ab13816f9808ed73307dff1efe227/data':
        'community.wave.seqera.io/library/typer:0.25.1--25ea8a9ce34456a3' }"

    input:
    tuple val(meta), path(files_in, stageAs: 'to_combine/*', arity: '1..*'), val(ids)

    output:
    tuple val(meta), path("*.saltshaker_classify.html"), emit: classify_html

    script:
    def prefix = task.ext.prefix ?: "${meta.caseid}"
    def args = task.ext.args ?: ""
    """
    saltshaker_to_html.py \
        $args \
        --input ${files_in.join(' --input ')} \
        --sample ${ids.join(' --sample ')} \
        --output ${prefix}.saltshaker_classify.html
    """
}
