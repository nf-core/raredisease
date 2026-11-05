process GET_CHROM_SIZES {
    tag "$fai"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c262fc09eca59edb5a724080eeceb00fb06396f510aefb229c2d2c6897e63975/data' :
        'community.wave.seqera.io/library/coreutils:9.5--ae99c88a9b28c264' }"

    input:
    tuple val(meta), path(fai)

    output:
    path '*.sizes'     , emit: sizes
    tuple val("${task.process}"), val('coreutils'), val("9.5"), topic: versions, emit: versions_coreutils

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cut -f 1,2 $fai > ${fai}.sizes
    """

    stub:
    """
    touch ${fai}.sizes
    """

}
