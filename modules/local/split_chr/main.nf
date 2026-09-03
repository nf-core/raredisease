process SPLIT_CHR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--h5bf99c6_6':
        'biocontainers/bioawk:1.0--h5bf99c6_6' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("chr_split"), emit: output
    tuple val("${task.process}"), val('bioawk'), val("1.0"), emit: versions_bioawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p chr_split
    bioawk -v outdir="chr_split" 'BEGIN{RS=">"; FS="\\n"} NR>1 {fnme=outdir "/" \$1 ".fa"; print ">" \$0 > fnme; close(fnme); if (NR==25) exit;}' ${input}
    """

    stub:
    """
    mkdir -p chr_split
    touch chr_split/test1.fa
    touch chr_split/test2.fa
    """
}
