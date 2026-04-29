process CONCAT_SALTSHAKER {
    tag "$meta.id"
    label "process_low"

    input:
    tuple val(meta), path(txts)

    output:
    tuple val(meta), path("*.txt"), emit: txt

    script:
    """
    cat ${txts.join(' ')} > ${meta.id}.saltshaker_classify.txt
    """
}
