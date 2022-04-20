process SENTIEON_DRIVER {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    input:
    tuple val(meta), path(bam), path(bai), path(score), path(score_idx), path(recal_pre), path(recal_post)
    path fasta
    path fai
    path known_dbsnp
    path known_mills
    path known_indels

    output:
    tuple val(meta), path('*.bam')                  , emit: bam          , optional: true
    tuple val(meta), path('*.bai')                  , emit: bai          , optional: true
    tuple val(meta), path('*.cram')                 , emit: cram         , optional: true
    tuple val(meta), path('*.crai')                 , emit: crai         , optional: true
    tuple val(meta), path('*.vcf.gz')               , emit: vcf          , optional: true
    tuple val(meta), path('*.vcf.gz.tbi')           , emit: vcf_tbi      , optional: true
    tuple val(meta), path('*recal_data.table')      , emit: recal_pre    , optional: true
    tuple val(meta), path('*recal_data.table.post') , emit: recal_post   , optional: true
    tuple val(meta), path('*recal.csv')             , emit: recal_csv    , optional: true
    tuple val(meta), path('*mq_metrics.txt')        , emit: mq_metrics   , optional: true
    tuple val(meta), path('*qd_metrics.txt')        , emit: qd_metrics   , optional: true
    tuple val(meta), path('*gc_summary.txt')        , emit: gc_summary   , optional: true
    tuple val(meta), path('*gc_metrics.txt')        , emit: gc_metrics   , optional: true
    tuple val(meta), path('*aln_metrics.txt')       , emit: aln_metrics  , optional: true
    tuple val(meta), path('*is_metrics.txt')        , emit: is_metrics   , optional: true
    tuple val(meta), path('*dedup_metrics.txt')     , emit: metrics_dedup, optional: true
    tuple val(meta), path('*score.txt.gz')          , emit: score        , optional: true
    tuple val(meta), path('*score.txt.gz.tbi')      , emit: score_idx    , optional: true
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def input  = bam ? '-i ' + bam.sort().join(' -i ') : ''
    def ref    = fasta ? "-r $fasta" : ''
    def dbsnp  = known_dbsnp  ? "-k $known_dbsnp" : ''
    def mills  = known_mills  ? "-k $known_mills" : ''
    def indels = known_indels ? "-k $known_indels" : ''
    if (args.contains('--algo Haplotyper')) {
        if (known_dbsnp) {
            dbsnp = ''
            def args_list = args.split('--algo Haplotyper')
            args_list = [ args_list[0] ] + ["--algo Haplotyper -d $known_dbsnp"] + [ args_list[-1] ]
            args = args_list.join(' ')
        }
    }
    """
    source sentieon_init.sh SENTIEON_LICENSE_BASE64

    sentieon \\
        driver \\
        $ref \\
        -t $task.cpus \\
        $input \\
        $dbsnp \\
        $mills \\
        $indels \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix       = task.ext.prefix ?: "${meta.id}"
    """
    source sentieon_init.sh SENTIEON_LICENSE_BASE64

    touch ${prefix}.bam
    touch ${prefix}.bai
    touch ${prefix}.cram
    touch ${prefix}.crai
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.recal_data.table
    touch ${prefix}.recal_data.table.post
    touch ${prefix}.recal.csv
    touch ${prefix}.dedup_metrics.txt
    touch ${prefix}.score.txt.gz
    touch ${prefix}.score.txt.gz.tbi
    touch ${prefix}.mq_metrics.txt
    touch ${prefix}.qd_metrics.txt
    touch ${prefix}.gc_summary.txt
    touch ${prefix}.gc_metrics.txt
    touch ${prefix}.aln_metrics.txt
    touch ${prefix}.is_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
