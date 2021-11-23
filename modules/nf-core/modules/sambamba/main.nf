// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SAMBAMBA {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }


    conda (params.enable_conda ? 'bioconda::sambamba=0.8.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/sambamba:0.8.1--h41abebc_0" 
    } else {
        container "quay.io/biocontainers/sambamba:0.8.1--h41abebc_0"
    }

    input:
    tuple val(meta), path(bam), path(bai) 

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"            , emit: version
   // To be changed in the script to chrM 
    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    sambamba view -t $task.cpus $bam -f bam chr20 > ${prefix}_mito.bam  
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(sambamba --version 2>&1) | sed 's/^.*sambamba //; s/Using.*\$//')
    END_VERSIONS
    """
}

