# nf-core/raredisease: Usage

**We recommend reading this documentation on the nf-core website: [https://nf-co.re/raredisease/usage](https://nf-co.re/raredisease/usage)**

Table of contents:

- [nf-core/raredisease: Usage](#nf-core-raredisease--usage)
  - [Introduction](#introduction)
  - [Prerequisites](#prerequisites)
  - [Run nf-core/raredisease with test data](#run-nf-core-raredisease-with-test-data)
  - [Run nf-core/raredisease with your data](#run-nf-core-raredisease-with-your-data)
    - [Samplesheet](#samplesheet)
    - [Reference files and parameters](#reference-files-and-parameters)
      - [1. Alignment](#1-alignment)
      - [2. QC stats from the alignment files](#2-qc-stats-from-the-alignment-files)
      - [3. Repeat expansions](#3-repeat-expansions)
      - [4. Variant calling - SNV](#4-variant-calling---snv)
      - [5. Variant calling - Structural variants](#5-variant-calling---structural-variants)
      - [6. SNV annotation & Ranking](#6-snv-annotation---ranking)
      - [7. SV annotation & Ranking](#7-sv-annotation---ranking)
      - [8. Mitochondrial analysis](#8-mitochondrial-analysis)
    - [Run the pipeline](#run-the-pipeline)
      - [Direct input in cli](#direct-input-in-cli)
      - [Import from a config file (recommended)](#import-from-a-config-file--recommended-)
  - [Best practices](#best-practices)
  - [Troubleshooting](#troubleshooting)
    - [Resource errors](#resource-errors)
      - [For beginners](#for-beginners)
      - [Advanced option on process level](#advanced-option-on-process-level)
  - [Custom configuration](#custom-configuration)
    - [Updating containers (advanced users)](#updating-containers--advanced-users-)
    - [nf-core/configs](#nf-core-configs)
    - [Run sentieon](#run-sentieon)
    - [Azure Resource Requests](#azure-resource-requests)
    - [Running in the background](#running-in-the-background)
    - [Nextflow memory requirements](#nextflow-memory-requirements)

## Introduction

nf-core/raredisease is a bioinformatics best-practice analysis pipeline to call, annotate and score variants from WGS/WES of rare disease patients. The pipeline is built using Nextflow.

## Prerequisites

1. Install Nextflow (>=22.10.1) using the instructions [here.](https://nextflow.io/docs/latest/getstarted.html#installation)
2. Install one of the following technologies for full pipeline reproducibility: Docker, Singularity, Podman, Shifter or Charliecloud.
   > Almost all nf-core pipelines give you the option to use conda as well. However, some tools used in the raredisease pipeline do not have a conda package so we do not support conda at the moment.

## Run nf-core/raredisease with test data

Before running the pipeline with your data, we recommend running it with the test dataset available [here](https://github.com/nf-core/test-datasets/tree/raredisease).

> You do not need to download the data as the pipeline is configured to fetch that data automatically for you when you use the test profile.

Run the following command, where YOURPROFILE is the package manager you installed on your machine. For example, `-profile test,docker` or `-profile test,singularity`:

```
nextflow run nf-core/raredisease \
    -revision dev -profile test,<YOURPROFILE> \
    --outdir <OUTDIR>
```

> Check [nf-core/configs](https://github.com/nf-core/configs/tree/master/conf) to see if a custom config file to run nf-core pipelines already exists for your institute. If so, you can simply use `-profile test,<institute>` in your command. This enables the appropriate package manager and sets the appropriate execution settings for your machine.
> NB: The order of profiles is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

The above command downloads the pipeline from GitHub, caches it, and tests it on the test dataset.

> When you run the command again, it will fetch the pipeline from cache even if a more recent version of the pipeline is available. To make sure that you're running the latest version of the pipeline, update the cached version of the pipeline by including `-latest` in the command.

Test profile runs the pipeline with a case containing three samples, but if you would like to test the pipeline with one sample, use `-profile test_one_sample,<YOURPROFILE>`.

Running the command creates the following files in your working directory:

```
work                # Directory containing the Nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other Nextflow hidden files, like history of pipeline logs.
```

## Run nf-core/raredisease with your data

Running the pipeline involves three steps:

1. Prepare a samplesheet
2. Gather all required references
3. Supply samplesheet and references, and run the command

#### Samplesheet

A samplesheet is used to pass the information about the sample(s), such as the path to the FASTQ files and other meta data (gender, phenotype, etc.,) to the pipeline in csv format.

nf-core/raredisease will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The pedigree information in the samplesheet (sex/gender and phenotype) should be provided as they would be for a [ped file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format) (i.e. 1 for male, 2 for female, other for unknown).

| Fields        | Description                                                                                                                                                                            |
| ------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`      | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `lane`        | Used to generate separate channels during the alignment step.                                                                                                                          |
| `fastq_1`     | Absolute path to FASTQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                         |
| `fastq_2`     | Absolute path to FASTQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                         |
| `gender`      | Sex (1=male; 2=female; other=unknown).                                                                                                                                                 |
| `phenotype`   | Affected status of patient (0 = missing; 1=unaffected; 2=affected).                                                                                                                    |
| `paternal_id` | Sample ID of the father, can be blank if the father isn't part of the analysis or for samples other than the proband.                                                                  |
| `maternal_id` | Sample ID of the mother, can be blank if the mother isn't part of the analysis or for samples other than the proband.                                                                  |
| `case_id`     | Case ID, for the analysis used when generating a family VCF.                                                                                                                           |

It is also possible to include multiple runs of the same sample in a samplesheet. For example, when you have re-sequenced the same sample more than once to increase sequencing depth. In that case, the `sample` identifiers in the samplesheet have to be the same. The pipeline will align the raw read/read-pairs independently before merging the alignments belonging to the same sample. Below is an example for a trio with the proband sequenced across two lanes:

| sample   | lane | fastq_1                          | fastq_2                          | gender | phenotype | paternal_id | maternal_id | case_id |
| -------- | ---- | -------------------------------- | -------------------------------- | ------ | --------- | ----------- | ----------- | ------- |
| AEG588A1 | 2    | AEG588A1_S1_L002_R1_001.fastq.gz | AEG588A1_S1_L002_R2_001.fastq.gz | 1      | 2         | AEG588A3    | AEG588A2    | fam_1   |
| AEG588A1 | 3    | AEG588A1_S1_L003_R1_001.fastq.gz | AEG588A1_S1_L003_R2_001.fastq.gz | 1      | 2         | AEG588A3    | AEG588A2    | fam_1   |
| AEG588A2 | 4    | AEG588A2_S1_L004_R1_001.fastq.gz | AEG588A2_S1_L004_R2_001.fastq.gz | 2      | 1         |             |             | fam_1   |
| AEG588A3 | 4    | AEG588A3_S1_L004_R1_001.fastq.gz | AEG588A3_S1_L004_R2_001.fastq.gz | 1      | 1         |             |             | fam_1   |

If you would like to see more examples of what a typical samplesheet looks like for a singleton and a trio, follow these links, [singleton](https://github.com/nf-core/test-datasets/blob/raredisease/testdata/samplesheet_single.csv) and [trio](https://github.com/nf-core/test-datasets/blob/raredisease/testdata/samplesheet_trio.csv).

#### Reference files and parameters

In nf-core/raredisease, references can be supplied using parameters listed [here](https://nf-co.re/raredisease/dev/parameters).

Note that the pipeline is modular in architecture. It offers you the flexibility to choose between different tools. For example, you can align with either bwamem2 or Sentieon BWA mem and call SNVs with either DeepVariant or Sentieon DNAscope. You also have the option to turn off sections of the pipeline if you do not want to run the. For example, snv annotation can be turned off by adding `--skip_snv_annotation` flag in the command line, or by setting it to true in a parameter file. This flexibility means that in any given analysis run, a combination of tools included in the pipeline will not be executed. So the pipeline is written in a way that can account for these differences while working with reference parameters. If a tool is not going to be executed during the course of a run, parameters used only by that tool need not be provided. For example, for SNV calling if you use DeepVariant as your variant caller, you need not provide the parameter `--ml_model`, which is only used by Sentieon DNAscope.

nf-core/raredisease consists of several tools used for various purposes. For convenience, we have grouped those tools under the following categories:

1. Alignment (bwamem2/Sentieon BWA mem)
2. QC stats from the alignment files
3. Repeat expansions (ExpansionsHunter & Stranger)
4. Variant calling - SNV (DeepVariant/Sentieon DNAscope)
5. Variant calling - Structural variants (SV) (Tiddit & Manta)
6. SNV annotation & ranking (rohcall, vcfanno, ensembl VEP, GENMOD)
7. SV annotation & ranking (SVDB query, ensembl VEP, GENMOD)
8. Mitochondrial analysis

> We have only listed the groups that require at least one input from the user. For example, the pipeline also runs SMNCopyNumberCaller, but it does not require any input other than the bam files passed by the pipeline. Hence, it is not mentioned in the list above. To know more about the tools used in the pipeline check the [README](../README.md).

The mandatory and optional parameters for each category are tabulated below.

> Alignment, QC stats, repeat expansions, SNV & SV variant calling are run by default. Hence, the mandatory parameters used by those features will always have to be provided to the pipeline.

##### 1. Alignment

| Mandatory           | Optional                    |
| ------------------- | --------------------------- |
| aligner<sup>1</sup> | fasta_fai<sup>2</sup>       |
| fasta               | bwamem2<sup>2</sup>         |
| platform            | known_dbsnp<sup>3</sup>     |
|                     | known_dbsnp_tbi<sup>3</sup> |

<sup>1</sup>Default value is bwamem2, but if you have a valid license for Sentieon, you have the option to use Sentieon as well.<br />
<sup>2</sup>fasta_fai and bwamem2, if not provided by the user, will be generated by the pipeline when necessary.<br />
<sup>3</sup>Used only by Sentieon.<br />

##### 2. QC stats from the alignment files

| Mandatory                                                    | Optional |
| ------------------------------------------------------------ | -------- |
| intervals_wgs<sup>1</sup>                                    |          |
| intervals_y<sup>1</sup>                                      |          |
| target_bed / (bait_intervals & target_intervals)<sup>2</sup> |          |

<sup>1</sup>These files are Picard's style interval list files, comprising the entire genome or only the chromosome Y. A version of these files for GRCh37 and for GRCh38 is supplied in the pipeline assets. These files are not necessary if you are using Sentieon.
<sup>2</sup> If a target_bed file is provided, bait_intervals and target_intervals can be generated by the pipeline.

##### 3. Repeat expansions

| Mandatory       | Optional |
| --------------- | -------- |
| variant_catalog |          |

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

##### 6. SNV annotation & Ranking

| Mandatory                     | Optional                       |
| ----------------------------- | ------------------------------ |
| genome<sup>1</sup>            | gnomad_af<sup>4</sup>          |
| vcfanno_resources<sup>2</sup> | reduced_penetrance<sup>5</sup> |
| vcfanno_toml<sup>3</sup>      | vcfanno_lua                    |
| vep_cache_version             | vep_filters<sup>6</sup>        |
| vep_cache                     | score_config_snv<sup>7</sup>   |
|                               | cadd_annotation<sup>8</sup>    |

<sup>1</sup>Genome version is used by VEP. You have the option to choose between GRCh37 and GRCh38.<br />
<sup>2</sup>Path to VCF files and their indices used by vcfanno. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/vcfanno_resources.txt).<br />
<sup>3</sup>Path to a vcfanno configuration file. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/vcfanno_config.toml).<br />
<sup>4</sup>GnomAD VCF file can be downloaded from [here] (https://gnomad.broadinstitute.org/downloads).<br />
<sup>5</sup>Used by GENMOD while modeling the variants. Contains a list of loci that show [reduced penetrance](https://medlineplus.gov/genetics/understanding/inheritance/penetranceexpressivity/) in people. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/reduced_penetrance.tsv).<br />
<sup>6</sup> This file contains a list of candidate genes (with [HGNC](https://www.genenames.org/) IDs) that is used to split the variants into canditate variants and research variants. Research variants contain all the variants, while candidate variants are a subset of research variants and are associated with candidate genes. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/hgnc.txt).<br />
<sup>7</sup>Used by GENMOD for ranking the variants. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/rank_model_snv.ini).<br />
<sup>8</sup>Path to a folder containing cadd annotations. Equivalent of the data/annotations/ folder described [here](https://github.com/kircherlab/CADD-scripts/#manual-installation), and it is used to calculate CADD scores for small indels. <br />

> NB: We use CADD only to annotate small indels. To annotate SNVs with precomputed CADD scores, pass the file containing CADD scores as a resource to vcfanno instead. Files containing the precomputed CADD scores for SNVs can be downloaded from [here](https://cadd.gs.washington.edu/download) (description: "All possible SNVs of GRCh3<7/8>/hg3<7/8>")

##### 7. SV annotation & Ranking

| Mandatory                  | Optional           |
| -------------------------- | ------------------ |
| genome                     | reduced_penetrance |
| svdb_query_dbs<sup>1</sup> | score_config_sv    |
| vep_cache_version          | vep_filters        |
| vep_cache                  |                    |

<sup>1</sup> A CSV file that describes the databases (VCFs) used by SVDB for annotating structural variants. Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/svdb_querydb_files.csv). Information about the column headers can be found [here](https://github.com/J35P312/SVDB#Query).

##### 8. Mitochondrial analysis

| Mandatory                      | Optional |
| ------------------------------ | -------- |
| genome                         |          |
| mt_backchain_shift<sup>1</sup> |          |
| mito_name                      |          |
| mt_fasta_shift                 |          |
| mt_intervals                   |          |
| mt_intervals_shift             |          |
| vcfanno_resources              |          |
| vcfanno_toml                   |          |
| vep_cache_version              |          |
| vep_cache                      |          |

<sup>1</sup>Can be generated by GATK's [ShiftFasta](https://gatk.broadinstitute.org/hc/en-us/articles/9570501436827-ShiftFasta-BETA-). Sample file [here](https://github.com/nf-core/test-datasets/blob/raredisease/reference/mt_shift8000.back_chain).

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

- **Restart a previous run:** Add `-resume` to your command when restarting a pipeline. Nextflow will use cached results from any pipeline steps where inputs are the same, and resume the run from where it terminated previously. For input to be considered the same, names and the files' contents must be identical. For more info about `-resume`, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html). You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

## Troubleshooting

#### Resource errors

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/raredisease/blob/dev/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/raredisease pipeline is failing after multiple re-submissions of the `DEEPVARIANT` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RAREDISEASE:RAREDISEASE:CALL_SNV:CALL_SNV_DEEPVARIANT:DEEPVARIANT (earlycasualcaiman)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RAREDISEASE:RAREDISEASE:CALL_SNV:CALL_SNV_DEEPVARIANT:DEEPVARIANT (earlycasualcaiman)'

Caused by:
    Process `NFCORE_RAREDISEASE:RAREDISEASE:CALL_SNV:CALL_SNV_DEEPVARIANT:DEEPVARIANT (earlycasualcaiman)` terminated with an error exit status (137)

Command executed:
    /opt/deepvariant/bin/run_deepvariant \
        --ref=reference.fasta \
        --reads=earlycasualcaiman_sorted_md.bam \
        --output_vcf=earlycasualcaiman_deepvar.vcf.gz \
        --output_gvcf=earlycasualcaiman_deepvar.g.vcf.gz \
        --model_type=WGS \
         \
        --num_shards=2

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed        /opt/deepvariant/bin/run_deepvariant --ref=reference.fasta     --reads=earlycasualcaiman_sorted_md.bam --output_vcf=earlycasualcaiman_deepvar.vcf.gz --output_gvcf=earlycasualcaiman_deepvar.g.vcf.gz --model_type=WGS --num_shards=2

Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

##### For beginners

nf-core/raredisease provides you with options to increase the amount of CPUs (`--max_cpus`), memory (`--max_memory`), and time (`--max_time`) for the whole pipeline. Based on the error above, you have to increase the amount of memory. Therefore you can go to the [parameter documentation of raredisease](https://nf-co.re/raredisease/dev/parameters) and click `show hidden params` button in the right panel to get the default value for `--max_memory`. Run the pipeline with updated max_memory value `--max_memory xxxGB` and `-resume` to skip all processes that were already calculated. If you cannot increase the resource of the complete pipeline, you can try to adapt the resource for a single process as mentioned below.

##### Advanced option on process level

To bypass this error you would need to find exactly which resources are set by the `DEEPVARIANT` process. The quickest way is to search for `process DEEPAVARIANT` in the [nf-core/raredisease Github repo](https://github.com/nf-core/raredisease/search?q=process+DEEPVARIANT).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/deepvariant/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_medium`](https://github.com/nf-core/raredisease/blob/7a0d47aca1d5b59771af2ce49c320249e379fc23/modules/nf-core/deepvariant/main.nf#L3).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default memory value for the `process_medium` label is set in the pipeline's [`base.config`](https://github.com/nf-core/raredisease/blob/7a0d47aca1d5b59771af2ce49c320249e379fc23/conf/base.config#L39-L43), which in this case is 36GB.
We can try and bypass the `DEEPVARIANT` process failure by creating a custom config file that increases the memory limit from 36GB to 72GB (NB: verify that the value set by --max_memory is above 72GB). The custom config below can then be provided to the pipeline via the `-c` parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RAREDISEASE:RAREDISEASE:CALL_SNV:CALL_SNV_DEEPVARIANT:DEEPVARIANT' {
        memory = 72.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_RAREDISEASE:RAREDISEASE:CALL_SNV:CALL_SNV_DEEPVARIANT:DEEPVARIANT` in the config file because this takes priority over the short name (`DEEPVARIANT`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

## Custom configuration

#### Updating containers (advanced users)

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

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
