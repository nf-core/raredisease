process PARSE_CONTAMINATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(contamination_table)

    output:
    tuple val(meta), path("*_contamination_mqc.tsv"), emit: mqc_table
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    
    import csv
    
    # Read GATK contamination table
    with open("${contamination_table}", 'r') as f:
        lines = f.readlines()
        # Skip header, get contamination value
        data_line = lines[1].strip().split('\\t')
        sample = data_line[0]
        contamination = float(data_line[1])
        contamination_pct = contamination * 100
    
    # Write MultiQC custom content file
    with open("${prefix}_contamination_mqc.tsv", 'w') as out:
        # Header with MultiQC configuration
        out.write("# id: 'gatk_contamination'\\n")
        out.write("# section_name: 'GATK Contamination'\\n")
        out.write("# description: 'Sample contamination estimates from GATK CalculateContamination'\\n")
        out.write("# plot_type: 'generalstats'\\n")
        out.write("# pconfig:\\n")
        out.write("#     contamination_pct:\\n")
        out.write("#         title: 'Contamination'\\n")
        out.write("#         description: 'Estimated sample contamination percentage'\\n")
        out.write("#         max: 10\\n")
        out.write("#         min: 0\\n")
        out.write("#         scale: 'RdYlGn-rev'\\n")
        out.write("#         suffix: '%'\\n")
        out.write("#         format: '{:,.2f}'\\n")
        # Data
        out.write("Sample\\tcontamination_pct\\n")
        out.write(f"${meta.id}\\t{contamination_pct:.4f}\\n")
    
    # Create versions file
    with open("versions.yml", 'w') as v:
        v.write('"${task.process}":\\n')
        v.write('    python: "3.11"\\n')
    """
}