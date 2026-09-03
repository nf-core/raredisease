process CREATE_PEDIGREE_FILE {
    tag "pedigree"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c262fc09eca59edb5a724080eeceb00fb06396f510aefb229c2d2c6897e63975/data' :
        'community.wave.seqera.io/library/coreutils:9.5--ae99c88a9b28c264' }"

    input:
    val(samples)

    output:
    path("*.ped"),       emit: ped
    tuple val("${task.process}"), val('coreutils'), val("9.5"), topic: versions, emit: versions_coreutils

    when:
    task.ext.when == null || task.ext.when

    script:
    def case_name = samples[0].case_id
    outfile_text = ['#family_id', 'sample_id', 'father', 'mother', 'sex', 'phenotype'].join('\\t')
    def samples_list = []
    samples.each { sample ->
        def sample_name =  sample.sample
        if (!samples_list.contains(sample_name)) {
            outfile_text += "\\n" + [sample.case_id, sample_name, sample.paternal, sample.maternal, sample.sex, sample.phenotype].join('\\t')
            samples_list.add(sample_name)
        }
    }
    """
    echo -e "$outfile_text" >${case_name}.ped
    """

    stub:
    def case_name = samples[0].case_id
    """
    touch ${case_name}.ped
    """
}
