process SVDB_QUERY {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::svdb=2.6.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svdb:2.6.1--py39h5371cbf_0':
        'quay.io/biocontainers/svdb:2.6.1--py39h5371cbf_0' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(in_occs), val(in_frqs), val(out_occs), val(out_frqs), path (vcf_dbs)
    val(placeholder)

    output:
    tuple val(meta), path("*_query.vcf")                                                               ,    emit: vcf
    tuple val(in_occs_rem), val(in_frqs_rem), val(out_occs_rem), val(out_frqs_rem), path (vcf_dbs_rem) ,    emit: vals
    path "versions.yml"                                                                                ,    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def in_occ        = "--in_occ "  + in_occs.head()
    in_occs_rem   = in_occs.tail()
    def in_frq        = "--in_frq "  + in_frqs.head()
    in_frqs_rem   = in_frqs.tail()
    def out_occ       = "--out_occ " + out_occs.head()
    out_occs_rem  = out_occs.tail()
    def out_frq       = "--out_frq " + out_frqs.head()
    out_frqs_rem  = out_frqs.tail()
    def vcf_db        = vcf_dbs.head()
    vcf_dbs_rem   = vcf_dbs.tail()
    """
    svdb \\
        --query \\
        $in_occ \\
        $in_frq \\
        $out_occ \\
        $out_frq \\
        $args \\
        --db $vcf_db \\
        --query_vcf $vcf \\
        --prefix ${prefix}_${task.index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_query.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
    END_VERSIONS
    """
}
