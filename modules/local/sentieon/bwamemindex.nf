process SENTIEON_BWAINDEX {
    tag "$fasta"
    label 'process_high'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    input:
    path fasta

    output:
    path "bwa_index"   , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sentieon_exe = params.sentieon_install_dir ? "${params.sentieon_install_dir}/sentieon" : 'sentieon'
    """
    source sentieon_init.sh SENTIEON_LICENSE_BASE64

    mkdir bwa_index

    $sentieon_exe bwa index \\
        $args \\
        -p bwa_index/${fasta.baseName} \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$($sentieon_exe driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$($sentieon_exe bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    def sentieon_exe = params.sentieon_install_dir ? "${params.sentieon_install_dir}/sentieon" : 'sentieon'
    """
    mkdir bwa_index

    touch bwa_index/${fasta}.ann
    touch bwa_index/${fasta}.pac
    touch bwa_index/${fasta}.amb
    touch bwa_index/${fasta}.bwt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$($sentieon_exe driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
        bwa: \$(echo \$($sentieon_exe bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
