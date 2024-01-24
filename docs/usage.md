# nf-core/raredisease: Usage

**We recommend reading this documentation on the nf-core website: [https://nf-co.re/raredisease/usage](https://nf-co.re/raredisease/usage)**

Table of contents:

- [nf-core/raredisease: Usage](#nf-coreraredisease-usage)
  - [Introduction](#introduction)
  - [Prerequisites](#prerequisites)
  - [Run nf-core/raredisease with test data](#run-nf-coreraredisease-with-test-data)
    - [Updating the pipeline](#updating-the-pipeline)
  - [Run nf-core/raredisease with your data](#run-nf-coreraredisease-with-your-data)
    - [Samplesheet](#samplesheet)
    - [Reference files and parameters](#reference-files-and-parameters)
      - [1. Alignment](#1-alignment)
      - [2. QC stats from the alignment files](#2-qc-stats-from-the-alignment-files)
      - [3. Repeat expansions](#3-repeat-expansions)
      - [4. Variant calling - SNV](#4-variant-calling---snv)
      - [5. Variant calling - Structural variants](#5-variant-calling---structural-variants)
      - [6. Copy number variant calling](#6-copy-number-variant-calling)
      - [7. SNV annotation \& Ranking](#7-snv-annotation--ranking)
      - [8. SV annotation \& Ranking](#8-sv-annotation--ranking)
      - [9. Mitochondrial annotation](#9-mitochondrial-annotation)
      - [10. Mobile element annoation](#10-mobile-element-annotation)
      - [11. Variant evaluation](#11-variant-evaluation)
    - [Run the pipeline](#run-the-pipeline)
      - [Direct input in CLI](#direct-input-in-cli)
      - [Import from a config file (recommended)](#import-from-a-config-file-recommended)
  - [Best practices](#best-practices)
  - [Core Nextflow arguments](#core-nextflow-arguments)
    - [`-profile`](#-profile)
    - [`-resume`](#-resume)
    - [`-c`](#-c)
  - [Custom configuration](#custom-configuration)
    - [Changing resources](#changing-resources)
    - [Custom Containers](#custom-containers)
    - [Custom Tool Arguments](#custom-tool-arguments)
      - [nf-core/configs](#nf-coreconfigs)
    - [Run Sentieon](#run-sentieon)
    - [Azure Resource Requests](#azure-resource-requests)
    - [Running in the background](#running-in-the-background)
    - [Nextflow memory requirements](#nextflow-memory-requirements)
    - [Running the pipeline without Internet access](#running-the-pipeline-without-internet-access)

## Introduction

nf-core/raredisease is a bioinformatics best-practice analysis pipeline to call, annotate and score variants from WGS/WES of rare disease patients. The pipeline is built using Nextflow.

## Prerequisites

1. Install Nextflow (>=22.10.1) using the instructions [here.](https://nextflow.io/docs/latest/getstarted.html#installation)
2. Install one of the following technologies for full pipeline reproducibility: Docker, Singularity, Podman, Shifter or Charliecloud.
   > Almost all nf-core pipelines give you the option to use conda as well. However, some tools used in the raredisease pipeline do not have a conda package so we do not support conda at the moment.

## Run nf-core/raredisease with test data

Before running the pipeline with your data, we recommend running it with the test dataset available [here](https://github.com/nf-core/test-datasets/tree/raredisease). You do not need to download the data as the pipeline is configured to fetch that data automatically for you when you use the test profile.

Run the following command, where YOURPROFILE is the package manager you installed on your machine. For example, `-profile test,docker` or `-profile test,singularity`:

```
nextflow run nf-core/raredisease \
    -revision dev -profile test,<YOURPROFILE> \
    --outdir <OUTDIR>
```

> Check [nf-core/configs](https://github.com/nf-core/configs/tree/master/conf) to see if a custom config file to run nf-core pipelines already exists for your institute. If so, you can simply use `-profile test,<institute>` in your command. This enables the appropriate package manager and sets the appropriate execution settings for your machine.
> NB: The order of profiles is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

Running the command creates the following files in your working directory:

```
work                # Directory containing the Nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other Nextflow hidden files, like history of pipeline logs.
```

Test profile runs the pipeline with a case containing three samples, but if you would like to test the pipeline with one sample, use `-profile test_one_sample,<YOURPROFILE>`.

> Note that the default cpu and memory configurations used in raredisease are written keeping the test profile (&dataset, which is tiny) in mind. You should override these values in configs to get it to work on larger datasets. Check the section `custom-configuration` below to know more about how to configure resources for your platform.

### Updating the pipeline

The above command downloads the pipeline from GitHub, caches it, and tests it on the test dataset. When you run the command again, it will fetch the pipeline from cache even if a more recent version of the pipeline is available. To make sure that you're running the latest version of the pipeline, update the cached version of the pipeline by including `-latest` in the command.

## Run nf-core/raredisease with your data

Running the pipeline involves three steps:

1. Prepare a samplesheet
2. Gather all required references
3. Supply samplesheet and references, and run the command

#### Samplesheet

A samplesheet is used to pass the information about the sample(s), such as the path to the FASTQ files and other meta data (sex, phenotype, etc.,) to the pipeline in csv format.

nf-core/raredisease will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The pedigree information in the samplesheet (sex and phenotype) should be provided as they would be for a [ped file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format) (i.e. 1 for male, 2 for female, other for unknown).

| Fields        | Description                                                                                                                                                                            |
| ------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`      | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `lane`        | Used to generate separate channels during the alignment step.                                                                                                                          |
| `fastq_1`     | Absolute path to FASTQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                         |
| `fastq_2`     | Absolute path to FASTQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                         |
| `sex`         | Sex (1=male; 2=female; other=unknown).                                                                                                                                                 |
| `phenotype`   | Affected status of patient (0 = missing; 1=unaffected; 2=affected).                                                                                                                    |
| `paternal_id` | Sample ID of the father, can be blank if the father isn't part of the analysis or for samples other than the proband.                                                                  |
| `maternal_id` | Sample ID of the mother, can be blank if the mother isn't part of the analysis or for samples other than the proband.                                                                  |
| `case_id`     | Case ID, for the analysis used when generating a family VCF.                                                                                                                           |

It is also possible to include multiple runs of the same sample in a samplesheet. For example, when you have re-sequenced the same sample more than once to increase sequencing depth. In that case, the `sample` identifiers in the samplesheet have to be the same. The pipeline will align the raw read/read-pairs independently before merging the alignments belonging to the same sample. Below is an example for a trio with the proband sequenced across two lanes:

| sample   | lane | fastq_1                          | fastq_2                          | sex | phenotype | paternal_id | maternal_id | case_id |
| -------- | ---- | -------------------------------- | -------------------------------- | --- | --------- | ----------- | ----------- | ------- |
| AEG588A1 | 2    | AEG588A1_S1_L002_R1_001.fastq.gz | AEG588A1_S1_L002_R2_001.fastq.gz | 1   | 2         | AEG588A3    | AEG588A2    | fam_1   |
| AEG588A1 | 3    | AEG588A1_S1_L003_R1_001.fastq.gz | AEG588A1_S1_L003_R2_001.fastq.gz | 1   | 2         | AEG588A3    | AEG588A2    | fam_1   |
| AEG588A2 | 4    | AEG588A2_S1_L004_R1_001.fastq.gz | AEG588A2_S1_L004_R2_001.fastq.gz | 2   | 1         |             |             | fam_1   |
| AEG588A3 | 4    | AEG588A3_S1_L004_R1_001.fastq.gz | AEG588A3_S1_L004_R2_001.fastq.gz | 1   | 1         |             |             | fam_1   |

If you would like to see more examples of what a typical samplesheet looks like for a singleton and a trio, follow these links, [singleton](https://github.com/nf-core/test-datasets/blob/raredisease/testdata/samplesheet_single.csv) and [trio](https://github.com/nf-core/test-datasets/blob/raredisease/testdata/samplesheet_trio.csv).

#### Reference files and parameters

In nf-core/raredisease, references can be supplied using parameters listed [here](https://nf-co.re/raredisease/dev/parameters).

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

Note that the pipeline is modular in architecture. It offers you the flexibility to choose between different tools. For example, you can align with bwamem2 or bwa or Sentieon BWA mem and call SNVs with either DeepVariant or Sentieon DNAscope. You also have the option to turn off sections of the pipeline if you do not want to run the. For example, snv annotation can be turned off by adding `--skip_snv_annotation` flag in the command line, or by setting it to true in a parameter file. This flexibility means that in any given analysis run, a combination of tools included in the pipeline will not be executed. So the pipeline is written in a way that can account for these differences while working with reference parameters. If a tool is not going to be executed during the course of a run, parameters used only by that tool need not be provided. For example, for SNV calling if you use DeepVariant as your variant caller, you need not provide the parameter `--ml_model`, which is only used by Sentieon DNAscope.

nf-core/raredisease consists of several tools used for various purposes. For convenience, we have grouped those tools under the following categories:

1. Alignment (bwamem2/bwa/Sentieon BWA mem)
2. QC stats from the alignment files
3. Repeat expansions (ExpansionsHunter & Stranger)
4. Variant calling - SNV (DeepVariant/Sentieon DNAscope)
5. Variant calling - Structural variants (SV) (Tiddit & Manta)
6. Copy number variant calling (GATK's GermlineCNVCaller)
7. SNV annotation & ranking (rohcall, vcfanno, ensembl VEP, GENMOD)
8. SV annotation & ranking (SVDB query, ensembl VEP, GENMOD)
9. Mitochondrial annotation

> We have only listed the groups that require at least one input from the user. For example, the pipeline also runs SMNCopyNumberCaller, but it does not require any input other than the bam files passed by the pipeline. Hence, it is not mentioned in the list above. To know more about the tools used in the pipeline check the [README](../README.md).

The mandatory and optional parameters for each category are tabulated below.

> Alignment, QC stats, repeat expansions, SNV & SV variant calling are run by default. Hence, the mandatory parameters used by those features will always have to be provided to the pipeline.

##### 1. Alignment

| Mandatory                      | Optional                       |
| ------------------------------ | ------------------------------ |
| aligner<sup>1</sup>            | fasta_fai<sup>4</sup>          |
| fasta<sup>2</sup>              | bwamem2<sup>4</sup>            |
| platform                       | bwa<sup>4</sup>                |
| mito_name/mt_fasta<sup>3</sup> | known_dbsnp<sup>5</sup>        |
|                                | known_dbsnp_tbi<sup>5</sup>    |
|                                | min_trimmed_length<sup>6</sup> |

<sup>1</sup>Default value is bwamem2. Other alternatives are bwa and sentieon (requires valid Sentieon license ).<br />
<sup>2</sup>Analysis set reference genome in fasta format, first 25 contigs need to be chromosome 1-22, X, Y and the mitochondria.<br />
<sup>3</sup>f If mito_name is provided, mt_fasta can be generated by the pipeline.<br />
<sup>4</sup>fasta_fai, bwa and bwamem2, if not provided by the user, will be generated by the pipeline when necessary.<br />
<sup>5</sup>Used only by Sentieon.<br />
<sup>6</sup>Default value is 40. Used only by fastp.<br />

##### 2. QC stats from the alignment files

| Mandatory                                                    | Optional |
| ------------------------------------------------------------ | -------- |
| intervals_wgs<sup>1</sup>                                    |          |
| intervals_y<sup>1</sup>                                      |          |
| target_bed / (bait_intervals & target_intervals)<sup>2</sup> |          |

<sup>1</sup>These files are Picard's style interval list files, comprising the entire genome or only the chromosome Y. A version of these files for GRCh37 and for GRCh38 is supplied in the pipeline assets. These files are not necessary if you are using Sentieon.<br />
<sup>2</sup> If a target_bed file is provided, bait_intervals and target_intervals can be generated by the pipeline.<br />

##### 3. Repeat expansions

| Mandatory                   | Optional |
| --------------------------- | -------- |
| variant_catalog<sup>1</sup> |          |

<sup>1</sup> We reccomend using the catalogs found [here](https://github.com/Clinical-Genomics/reference-files/tree/master/rare-disease/disease_loci/ExpansionHunter-v5.0.0). These catalogs have been extended from the illumina ones to include information on pathogenicity, which is neccesarry for the workflow.

##### 4. Variant calling - SNV

| Mandatory                  | Optional                    |
| -------------------------- | --------------------------- |
| variant_caller<sup>1</sup> | known_dbsnp<sup>2</sup>     |
| ml_model<sup>2</sup>       | known_dbsnp_tbi<sup>2</sup> |
| analysis_type<sup>3</sup>  | call_interval<sup>2</sup>   |
|                            | known_dbsnp_tbi<sup>2</sup> |

<sup>1</sup>Default variant caller is DeepVariant, but you have the option to use Sentieon as well.<br />
<sup>2</sup>These parameters are only used by Sentieon.<br />
<sup>3</sup>Default is WGS, but you have the option to choose WES as well.<br />

##### 5. Variant calling - Structural variants

| Mandatory | Optional   |
| --------- | ---------- |
|           | target_bed |
|           | bwa        |

##### 6. Copy number variant calling

| Mandatory                      | Optional                        |
| ------------------------------ | ------------------------------- |
| ploidy_model<sup>1</sup>       | readcount_intervals<sup>3</sup> |
| gcnvcaller_model<sup>1,2</sup> |                                 |

<sup>1</sup> Output from steps 3 & 4 of GATK's CNV calling pipeline run in cohort mode as described [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variants).<br />
<sup>2</sup> Sample file can be found [here](https://raw.githubusercontent.com/nf-core/test-datasets/raredisease/reference/gcnvmodels.tsv) (Note the header 'models' in the sample file).<br />
<sup>3</sup> Output from step 1 of GATK's CNV calling pipeline as described [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variants).<br />

##### 7. SNV annotation & Ranking

| Mandatory                     | Optional                       |
| ----------------------------- | ------------------------------ |
| genome<sup>1</sup>            | reduced_penetrance<sup>7</sup> |
| vcfanno_resources<sup>2</sup> | vcfanno_lua                    |
| vcfanno_toml<sup>3</sup>      | vep_filters<sup>8</sup>        |
| vep_cache_version             | cadd_resources<sup>9</sup>     |
| vep_cache<sup>4</sup>         | vep_plugin_files<sup>10</sup>  |
| gnomad_af<sup>5</sup>         |                                |
| score_config_snv<sup>6</sup>  |                                |

<sup>1</sup>Genome version is used by VEP. You have the option to choose between GRCh37 and GRCh38.<br />
<sup>2</sup>Path to VCF files and their indices used by vcfanno. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/vcfanno_resources.txt).<br />
<sup>3</sup>Path to a vcfanno configuration file. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/vcfanno_config.toml).<br />
<sup>4</sup> VEP caches can be downloaded [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache).
VEP plugins may be installed in the cache directory, and the plugin pLI is mandatory to install. To supply files required by VEP plugins, use `vep_plugin_files` parameter.
See example cache [here](https://raw.githubusercontent.com/nf-core/test-datasets/raredisease/reference/vep_cache_and_plugins.tar.gz).<br />
<sup>5</sup> GnomAD VCF files can be downloaded from [here](https://gnomad.broadinstitute.org/downloads). The option `gnomad_af` expects a tab-delimited file with
no header and the following columns: `CHROM POS REF_ALLELE ALT_ALLELE AF`. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/gnomad_reformated.tab.gz).<br />
<sup>6</sup>Used by GENMOD for ranking the variants. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/rank_model_snv.ini).<br />
<sup>7</sup>Used by GENMOD while modeling the variants. Contains a list of loci that show [reduced penetrance](https://medlineplus.gov/genetics/understanding/inheritance/penetranceexpressivity/) in people. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/reduced_penetrance.tsv).<br />
<sup>8</sup> This file contains a list of candidate genes (with [HGNC](https://www.genenames.org/) IDs) that is used to split the variants into canditate variants and research variants. Research variants contain all the variants, while candidate variants are a subset of research variants and are associated with candidate genes. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/hgnc.txt). Not required if --skip_vep_filter is set to true.<br />
<sup>9</sup>Path to a folder containing cadd annotations. Equivalent of the data/annotations/ folder described [here](https://github.com/kircherlab/CADD-scripts/#manual-installation), and it is used to calculate CADD scores for small indels. <br />
<sup>10</sup>A CSV file that describes the files used by VEP's named and custom plugins. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/vep_files.csv). <br />

> NB: We use CADD only to annotate small indels. To annotate SNVs with precomputed CADD scores, pass the file containing CADD scores as a resource to vcfanno instead. Files containing the precomputed CADD scores for SNVs can be downloaded from [here](https://cadd.gs.washington.edu/download) (description: "All possible SNVs of GRCh3<7/8>/hg3<7/8>")

##### 8. SV annotation & Ranking

| Mandatory                                      | Optional           |
| ---------------------------------------------- | ------------------ |
| genome                                         | reduced_penetrance |
| svdb_query_dbs/svdb_query_bedpedbs<sup>1</sup> |                    |
| vep_cache_version                              | vep_filters        |
| vep_cache                                      | vep_plugin_files   |
| score_config_sv                                |                    |

<sup>1</sup> A CSV file that describes the databases (VCFs or BEDPEs) used by SVDB for annotating structural variants. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/svdb_querydb_files.csv). Information about the column headers can be found [here](https://github.com/J35P312/SVDB#Query).

##### 9. Mitochondrial annotation

| Mandatory         | Optional         |
| ----------------- | ---------------- |
| genome            | vep_filters      |
| mito_name         | vep_plugin_files |
| vcfanno_resources |                  |
| vcfanno_toml      |                  |
| vep_cache_version |                  |
| vep_cache         |                  |
| score_config_mt   |                  |

##### 10. Mobile element annotation

| Mandatory                                   | Optional    |
| ------------------------------------------- | ----------- |
| genome                                      | vep_filters |
| mobile_element_svdb_annotations<sup>1</sup> |             |
| vep_cache_version                           |             |
| vep_cache                                   |             |

<sup>1</sup> A CSV file that describes the databases (VCFs) used by SVDB for annotating mobile elements with allele frequencies. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/svdb_querydb_files.csv).

##### 11. Variant evaluation

| Mandatory                  | Optional |
| -------------------------- | -------- |
| run_rtgvcfeval<sup>1</sup> | sdf      |
| rtg_truthvcfs<sup>2</sup>  |          |

<sup>1</sup> This parameter is set to false by default, set it to true if if you'd like to run the evaluation subworkflow
<sup>2</sup> A CSV file that describes the truth VCF files used by RTG Tools' vcfeval for evaluating SNVs. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/rtg_example.csv). The file contains four columns `samplename,vcf,bedregions,evaluationregions` where samplename is the user assigned samplename in the input samplesheet, vcf is the path to the truth vcf file, bedregions and evaluationregions are the path to the bed files that are supposed to be passed through --bed_regions and --evaluation_regions options of vcfeval.

#### Run the pipeline

You can directly supply the parameters in the command line (CLI) or use a config file from which the pipeline can import the parameters.

##### Direct input in CLI

All of the pipeline parameters listed [here](https://nf-co.re/raredisease/dev/parameters) can be passed using the CLI. For example, you can provide the samplesheet, reference FASTA, and turn off annotations for SNVs and SVs by running,

```
nextflow run nf-core/raredisease \
    -revision dev \
    -profile test,<YOURPROFILE> \
    --input samplesheet.csv \
    --fasta reference.fasta \
    --skip_snv_annotation \
    --skip_sv_annotation \
    --outdir <OUTDIR>
```

##### Import from a config file (recommended)

To input or change the default parameters, we recommend using a config file instead. Create a config file that contains all the parameters and supply that file as shown below.

```
nextflow run nf-core/raredisease \
    -revision dev \
    -profile test,<YOURPROFILE> \
    -c <YOURCONFIG> \
    --outdir <OUTDIR>
```

A sample config file can be found [here](https://github.com/nf-core/raredisease/blob/dev/conf/test.config).

## Best practices

- **Singularity cache:** If you are using singularity, use the nf-core download command to download images first, before running the pipeline. Define [NXF_SINGULARITY_CACHEDIR](https://nextflow.io/docs/latest/config.html?highlight=nxf_singularity_cachedir#environment-variables) or singularity.cacheDir Nextflow options to store and re-use the images from a central location for future pipeline runs.

- **Save references:** While the pipeline can generate indices for the FASTA and some VCF files when they are not provided, it can be benficial to supply them to the pipeline as it saves computing resources. You can use the `--save_reference` parameter to publish those files, and then update your config file with the paths to these files for your subsequent runs.

- **Update pipeline:** If you would like to update the pipeline code stored in cache, in addition to running the command with `-latest` option, you can also run the command below:

```bash
nextflow pull nf-core/raredisease
```

- **Reproducibility:** Specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since. First, go to the [nf-core/raredisease releases page](https://github.com/nf-core/raredisease/releases) and find the latest pipeline version - numeric only (e.g. `1.3.1`). Then specify this when running the pipeline with `-r`, for example, `-r 1.3.1`. You can switch to another version by changing the number after the `-r` flag. The version number will be logged in reports when you run the pipeline. For example, you can view the version number at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

- **Restart a previous run:** Add `-resume` to your command when restarting a pipeline. Nextflow will use cached results from any pipeline steps where inputs are the same, and resume the run from where it terminated previously. For input to be considered the same, names and the files' contents must be identical. For more info about `-resume`, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html). You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

- **Reusing parameters:** If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> ⚠️ Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
> The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/raredisease -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
input: 'data'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

{%- if nf_core_configs %}

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).
{% else %}
{% endif %}
Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Changing resources

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

#### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Run Sentieon

To use Sentieon you have to:

1. Ensure that Sentieon executable is in path.
2. If necessary, store license information as a secret in Nextflow's local store.
3. Set environmental variable `NXF_ENABLE_SECRETS` to an appropriate value.

To elaborate more on item #2 in the list above, if you are running nf-core/raredisease with Sentieon on AWS or on a local machine and you do not want other users to know your license details, we recommend that you use [Nextflow's secrets feature](https://www.nextflow.io/docs/latest/secrets.html) to pass the that information. To do that run the command below after replacing LICENSE with the value corresponding to your license server (for example, 1.1.1.1:4000)

```
nextflow secrets set SENTIEON_LICENSE_BASE64 <LICENSE>
```

If you are using Nextflow secrets, you have to set the environment variable `NXF_ENABLE_SECRETS` to true. This will see to it that the pipeline can retrieve the secret from Nextflow's secrets store during the pipeline execution. Keep in mind that versions of Nextflow Version 22.09.2-edge and onwards have NXF_ENABLE_SECRETS to true by default. If you are not using secrets set `NXF_ENABLE_SECRETS` to false, but make sure that the environment variable [`SENTIEON_LICENSE`](`NXF_ENABLE_SECRETS`) is set to reflect the value of your license server on your machine.

### Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run Nextflow within a cluster job submitted to your job scheduler (from where it submits more jobs).

### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

### Running the pipeline without Internet access

The pipeline and container images can be downloaded using [nf-core tools](https://nf-co.re/docs/usage/offline). For running offline, you of course have to make all the reference data available locally, and specify `--fasta`, etc., see [above](#reference-files-and-parameters).

Contrary to the paragraph about [Nextflow](https://nf-co.re/docs/usage/offline#nextflow) on the page linked above, it is not possible to use the "-all" packaged version of Nextflow for this pipeline. The online version of Nextflow is necessary to support the necessary nextflow plugins. Download instead the file called just `nextflow`. Nextflow will download its dependencies when it is run. Additionally, you need to download the nf-validation plugin explicitly:

```
./nextflow plugin install nf-validation
```

Now you can transfer the `nextflow` binary as well as its directory `$HOME/.nextflow` to the system without Internet access, and use it there. It is necessary to use an explicit version of `nf-validation` offline, or Nextflow will check for the most recent version online. Find the version of nf-validation you downloaded in `$HOME/.nextflow/plugins`, then specify this version for `nf-validation` in your configuration file:

```
plugins {
        // Set the plugin version explicitly, otherwise nextflow will look for the newest version online.
        id 'nf-validation@0.3.1'
}
```

This should go in your Nextflow confgiguration file, specified with `-c <YOURCONFIG>` when running the pipeline.
