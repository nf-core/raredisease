# nf-core/raredisease: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

### Alignment

#### Mapping

##### Bwa-mem2

[Bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) used to map the reads to a reference genome. The aligned reads are coordinate sorted with samtools sort. These files are treated as intermediates and are not placed in the output folder by default.

##### Sentieon

[Sentieon's bwa mem](https://support.sentieon.com/manual/DNAseq_usage/dnaseq/#map-reads-to-reference) is the software accelerated version of the bwa-mem algorithm. It is used to efficiently perform the alignment using BWA. Aligned reads are then coordinate sorted using Sentieon's [sort](https://support.sentieon.com/manual/usages/general/#util-syntax) utility. It is not the default aligner, but it can be chosen over bwamem2 by setting `--aligner` option to sentieon.

#### Duplicate marking

##### Picard's MarkDuplicates

[Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)

##### Sentieon dedup

[Sentieon dedup](https://support.sentieon.com/manual/DNAseq_usage/dnaseq/#remove-or-mark-duplicates)

<details markdown="1">
<summary>Output files from Alignment</summary>

- `{outputdir}/alignment/`
  - `*.bam`: FastQC report containing quality metrics.
  - `*.bai`: Zip archive containing the FastQC report, tab-delimited data file and plot images.
  - `*.txt`: Text file containing the dedup metrics.
  </details>

### Quality control and reporting

#### Quality control

##### FastQC

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

##### Mosdepth

[Mosdepth](https://github.com/brentp/mosdepth)

##### Picardtools

[Picard's CollectMutipleMetrics, CollectHsMetrics, and CollectWgsMetrics](https://broadinstitute.github.io/picard/)

##### Qualimap

[Qualimap](http://qualimap.conesalab.org/)

##### Sention WgsMetricsAlgo

[Sentieon's WgsMetricsAlgo](https://support.sentieon.com/manual/usages/general/)

#### TIDDIT cov

[TIDDIT's cov](https://github.com/J35P312/)

#### Reporting

##### MultiQC

[MultiQC](http://multiqc.info/)

### Variant calling - SNV

#### DeepVariant

[DeepVariant](https://github.com/google/deepvariant)

#### Sentieon DNAscope

[Sentieon DNAscope](https://support.sentieon.com/manual/DNAscope_usage/dnascope/)

### Variant calling - SV

#### Manta

[Manta](https://github.com/Illumina/manta)

#### TIDDIT sv

[TIDDIT's sv](https://github.com/SciLifeLab/TIDDIT)

### Variant calling - repeat expansions

#### ExpansionsHunter

[ExpansionHunter](https://github.com/Illumina/ExpansionHunter)

#### stranger

[stranger](https://github.com/Clinical-Genomics/stranger)

### Annotation - SNV

#### bcftools roh

[bcftools roh](https://samtools.github.io/bcftools/bcftools.html#roh)

#### VCFanno

[VCFanno](https://github.com/brentp/vcfanno)

#### VEP

[VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)

### Annotation - SV

#### VCFanno

[VCFanno](https://github.com/brentp/vcfanno)

#### VEP

[VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)

### Mitochondrial analysis

#### Alignment and variant calling

[Alignment and variant calling - GATK Mitochondrial short variant discovery pipeline ](https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels-)

#### Annotation:

##### Haplogrep2

[Haplogrep2](https://github.com/seppinho/haplogrep-cmd)

##### VCFanno

[VCFanno](https://github.com/brentp/vcfanno)

##### VEP

[VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)

### Rank variants and filtering

#### Genmod

[Genmod](https://github.com/Clinical-Genomics/genmod)

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
