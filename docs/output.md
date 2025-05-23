# nf-core/raredisease: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [nf-core/raredisease: Output](#nf-coreraredisease-output)
  - [Introduction](#introduction)
  - [Pipeline overview](#pipeline-overview)
    - [Alignment](#alignment)
      - [Mapping](#mapping)
        - [Bwa-mem2](#bwa-mem2)
        - [BWA](#bwa)
        - [BWA-MEME](#bwa-meme)
        - [Sentieon bwa mem](#sentieon-bwa-mem)
      - [Duplicate marking](#duplicate-marking)
        - [Picard's MarkDuplicates](#picards-markduplicates)
        - [Sentieon Dedup](#sentieon-dedup)
      - [Subsample mitochondrial alignments](#subsample-mitochondrial-alignments)
    - [Quality control and reporting](#quality-control-and-reporting)
      - [Quality control](#quality-control)
        - [FastQC](#fastqc)
        - [Mosdepth](#mosdepth)
        - [Picard tools](#picard-tools)
        - [Qualimap](#qualimap)
        - [Chromograph coverage](#chromograph-coverage)
        - [Sention WgsMetricsAlgo](#sention-wgsmetricsalgo)
        - [TIDDIT's cov and UCSC WigToBigWig](#tiddits-cov-and-ucsc-wigtobigwig)
        - [VerifyBamID2](#verifybamid2)
      - [Reporting](#reporting)
        - [MultiQC](#multiqc)
    - [Variant calling - SNV](#variant-calling---snv)
      - [DeepVariant](#deepvariant)
      - [Sentieon DNAscope](#sentieon-dnascope)
    - [Variant calling - SV](#variant-calling---sv)
      - [Manta](#manta)
      - [TIDDIT sv](#tiddit-sv)
      - [GATK GermlineCNVCaller - CNV calling](#gatk-germlinecnvcaller---cnv-calling)
      - [CNVnator - CNV calling](#cnvnator---cnv-calling)
      - [SVDB merge](#svdb-merge)
    - [Variant calling - repeat expansions](#variant-calling---repeat-expansions)
      - [Expansion Hunter](#expansion-hunter)
      - [Stranger](#stranger)
    - [Annotation - SNV](#annotation---snv)
      - [bcftools roh](#bcftools-roh)
      - [vcfanno](#vcfanno)
      - [CADD](#cadd)
      - [VEP](#vep)
      - [UPD](#upd)
      - [Chromograph](#chromograph)
      - [Rhocall viz](#rhocall-viz)
    - [Annotation - SV](#annotation---sv)
      - [SVDB query](#svdb-query)
      - [VEP](#vep-1)
    - [Mitochondrial analysis](#mitochondrial-analysis)
      - [Haplocheck](#haplocheck)
      - [Alignment and variant calling](#alignment-and-variant-calling)
        - [MT deletion script](#mt-deletion-script)
        - [eKLIPse](#eklipse)
      - [Annotation](#annotation)
        - [HaploGrep2](#haplogrep2)
        - [vcfanno](#vcfanno-1)
        - [CADD](#cadd-1)
        - [Hmtnote](#hmtnote)
        - [VEP](#vep-2)
    - [Filtering and ranking](#filtering-and-ranking)
      - [Filter_vep](#filter_vep)
      - [GENMOD](#genmod)
    - [Mobile element analysis](#mobile-element-analysis)
      - [Calling mobile elements](#calling-mobile-elements)
      - [Annotating mobile elements](#annotating-mobile-elements)
    - [Variant evaluation](#variant-evaluation)
    - [Gens](#gens)
    - [Pipeline information](#pipeline-information)

### Alignment

#### Mapping

##### Bwa-mem2

[Bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) used to map the reads to a reference genome. The aligned reads are coordinate sorted with samtools sort. These files are treated as intermediates and are not placed in the output folder by default.

##### BWA

[BWA](https://github.com/lh3/bwa) used to map the reads to a reference genome. The aligned reads are coordinate sorted with samtools sort. These files are treated as intermediates and are not placed in the output folder by default. It is not the default aligner, but it can be chosen by setting `--aligner` option to bwa.

##### BWA-MEME

[BWA-MEME](https://github.com/kaist-ina/BWA-MEME) used to map the reads to a reference genome. The aligned reads are coordinate sorted with samtools sort. These files are treated as intermediates and are not placed in the output folder by default. It is not the default aligner, but it can be chosen by setting `--aligner` option to bwameme.

##### Sentieon bwa mem

[Sentieon's bwa mem](https://support.sentieon.com/manual/DNAseq_usage/dnaseq/#map-reads-to-reference) is the software accelerated version of the bwa-mem algorithm. It is used to efficiently perform the alignment using BWA. Aligned reads are then coordinate sorted using Sentieon's [sort](https://support.sentieon.com/manual/usages/general/#util-syntax) utility. These files are treated as intermediates and are not placed in the output folder by default. It is not the default aligner, but it can be chosen by setting `--aligner` option to "sentieon".

#### Duplicate marking

##### Picard's MarkDuplicates

[Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) is used for marking PCR duplicates that can occur during library amplification. This is essential as the presence of such duplicates results in false inflated coverages, which in turn can lead to overly-confident genotyping calls during variant calling. Only reads aligned by Bwa-mem2, bwameme and bwa are processed by this tool. By default, alignment files are published in bam format. If you would like to store cram files instead, set `--save_mapped_as_cram` to true.

<details markdown="1">
<summary>Output files from Alignment</summary>

- `{outputdir}/alignment/`
  - `*.bam|*.cram`: Alignment file in bam/cram format.
  - `*.bai|*.crai`: Index of the corresponding bam/cram file.
  - `*.txt`: Text file containing the dedup metrics.
  </details>

##### Sentieon Dedup

[Sentieon Dedup](https://support.sentieon.com/manual/DNAseq_usage/dnaseq/#remove-or-mark-duplicates) is the algorithm used by Sentieon's driver to remove duplicate reads. Only reads aligned by Sentieon's implementation of bwa are processed by this algorithm. By default, alignment files are published in bam format. If you would like to store cram files instead, set `--save_mapped_as_cram` to true.

<details markdown="1">
<summary>Output files from Alignment</summary>

- `{outputdir}/alignment/`
  - `*.bam|*.cram`: Alignment file in bam/cram format.
  - `*.bai|*.crai`: Index of the corresponding bam/cram file.
  - `*.metrics`: Text file containing the dedup metrics.
  </details>

#### Subsample mitochondrial alignments

[Samtools view](https://www.htslib.org/doc/samtools-view.html) is used by the pipeline to subsample mitochondrial alignments to a user specified coverage. The file is mainly intended to be used for visualization of MT alignments in IGV. The non-subsampled bam file is used for variant calling and other downstream analysis steps.

<details markdown="1">
<summary>Output files from Alignment</summary>

- `{outputdir}/alignment/`
  - `<sampleid>_mt_subsample.bam`: Alignment file in bam format.
  - `<sampleid>_mt_subsample.bam.bai`: Index of the corresponding bam file.
  </details>

### Quality control and reporting

#### Quality control

##### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

<details markdown="1">
<summary>Output files</summary>

- `{outputdir}/fastqc/{sampleid}_T*/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

:::note
The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.
:::

##### Mosdepth

[Mosdepth](https://github.com/brentp/mosdepth) is used to report quality control metrics such as coverage, and GC content from alignment files. The global distribution file, generated by this program is passed to MultiQC for generating the following plots,

- Cumulative coverage distribution
- Coverage distribution
- Average coverage per contig

<details markdown="1">
<summary>Output files</summary>

- `{outputdir}/qc_bam/`
  - `<sampleid>_mosdepth.global.dist.txt`: This file contains a cumulative distribution indicating the proportion of total bases that were covered for at least a given coverage value across each chromosome and the whole genome.
  - `<sampleid>_mosdepth.summary.txt`: This file contains summary statistics, such as mean, minimum and maximum coverage per genomic contig.

</details>

##### Picard tools

[Picard's CollectMutipleMetrics, CollectHsMetrics, and CollectWgsMetrics](https://broadinstitute.github.io/picard/) We use Picardtools' CollectWgsMetrics and CollectHsMetrics utilities to calculate metrics about coverage and performance of WGS & WES experiments. In addition to those metrics, we use CollectMultipleMetrics to gather information about alignment summary, insert size, GC content etc., The metrics generated by these three utilites are passed along to MultiQC to generate several plots as well.

<details markdown="1">
<summary>Output files</summary>

- `{outputdir}/qc_bam/<sampleid>_qualimap/`
  - `<sampleid>_hsmetrics.CollectHsMetrics.coverage_metrics`:
  - `<sampleid>_multiplemetrics.CollectMultipleMetrics.alignment_summary_metrics`:
  - `<sampleid>_multiplemetrics.CollectMultipleMetrics.base_distribution_by_cycle_metrics`:
  - `<sampleid>_multiplemetrics.CollectMultipleMetrics.insert_size_metrics`:
  - `<sampleid>_multiplemetrics.CollectMultipleMetrics.quality_by_cycle_metrics`:
  - `<sampleid>_multiplemetrics.CollectMultipleMetrics.quality_distribution_metrics`:
  - `<sampleid>_wgsmetrics.CollectWgsMetrics.coverage_metrics`:
  - `<sampleid>_wgsmetrics_y.CollectWgsMetrics.coverage_metrics`:
  </details>

##### Qualimap

[Qualimap](http://qualimap.conesalab.org/) also allows you to assess the alignment coverage. Qualimap results are used by MultiQC to generate the following plots.

- Coverage histogram
- Cumulative genome coverage
- Insert size histogram
- GC content distribution

<details markdown="1">
<summary>Output files</summary>

- `{outputdir}/qc_bam/<sampleid>_qualimap/` this directory includes a qualimap report and associated raw statistic files. You can open the .html file in your internet browser to see the in-depth report.
  </details>

##### Chromograph coverage

[Chromograph](https://github.com/Clinical-Genomics/chromograph) is a python package to create PNG images from genetics data such as BED and WIG files.

<details markdown="1">
<summary>Output files</summary>

- `{outputdir}/qc_bam/<sampleid>_chromographcov/*.png` plots showing coverage across chromosomes for each chromosome.
  </details>

##### Sention WgsMetricsAlgo

[Sentieon's WgsMetricsAlgo](https://support.sentieon.com/manual/usages/general/) is the Sentieon's equivalent of Picard's CollectWgsMetrics.

<details markdown="1">
<summary>Output files</summary>

- `{outputdir}/qc_bam/`
  - `<sampleid>_wgsmetrics.txt`:
  </details>

##### TIDDIT's cov and UCSC WigToBigWig

[TIDDIT's cov](https://github.com/J35P312/) is used to analyse the read depth of a bam file and generates a coverage report in wig format. This file is later passed to [UCSC WigToBigWig](https://genome.ucsc.edu/goldenPath/help/bigWig.html) to convert the file into a bigwig.

<details markdown="1">
<summary>Output files</summary>

- `{outputdir}/qc_bam/`
  - `<sampleid>_tidditcov.wig`:
  - `<sampleid>.bw`:
    </details>

##### VerifyBamID2

[VerifyBamID2](https://github.com/Griffan/VerifyBamID) is used to analyse a bam file and generates a contamination report. The pipeline will only generate the following files when the parameters `verifybamid_svd_bed`, `verifybamid_svd_mu`, and `verifybamid_svd_ud` are provided.

<details markdown="1">
<summary>Output files</summary>

- `{outputdir}/qc_bam/`
  - `<sampleid>.selfSM`:
  - `<sampleid>.Ancestry`:
    </details>

#### Reporting

##### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

### Variant calling - SNV

#### DeepVariant

[DeepVariant](https://github.com/google/deepvariant) is a deep learning-based variant caller that takes aligned reads, produces pileup image tensors from them, classifies each tensor using a convolutional neural network and finally reports the results in a standard VCF or gVCF file. Variant calls generated by DeepVariant are joint genotyped using [GLnexus](https://github.com/dnanexus-rnd/GLnexus), and then normalized using bcftools norm. Only the normalized vcfs are placed in the output folder by default.

> **NB**: In case you are running the separate mitochondrial analysis, mitochondrial calls are filtered from the normalized vcfs before they are published using GATK SelectVariants.

<details markdown="1">
<summary>Output files</summary>

- `call_snv/genome`
  - `<case_id>_snv.vcf.gz`: normalized vcf file containing no MT variants.
  - `<case_id>_snv.vcf.gz.tbi`: index of the vcf file containing no MT variants.

</details>

#### Sentieon DNAscope

The pipeline performs variant calling using [Sentieon DNAscope](https://support.sentieon.com/manual/DNAscope_usage/dnascope/) with a machine learning model. This approach identifies the candidate sites with a higher accuracy, and calculates genotypes for each sample at that site. These files are treated as intermediates and are not placed in the output folder by default. DNAscope is not run by default. To use DNAscope instead of DeepVariant, set `--variant_caller` to sentieon.

<details markdown="1">
<summary>Output files</summary>

- `call_snv/genome`
  - `<case_id>_snv.vcf.gz`: normalized vcf file containing no MT variants.
  - `<case_id>_snv.vcf.gz.tbi`: index of the vcf file containing no MT variants.

</details>

### Variant calling - SV

#### Manta

[Manta](https://github.com/Illumina/manta) calls structural variants (SVs) and indels from mapped paired-end sequencing reads. It combines paired and split-read evidence during SV discovery and scoring to improve accuracy, but does not require split-reads or successful breakpoint assemblies to report a variant in cases where there is strong evidence otherwise. Output vcf files are treated as intermediates and are not placed in the output folder.

#### TIDDIT sv

[TIDDIT's sv](https://github.com/SciLifeLab/TIDDIT) is used to identify chromosomal rearrangements using sequencing data. TIDDIT identifies intra and inter-chromosomal translocations, deletions, tandem-duplications and inversions, using supplementary alignments as well as discordant pairs. TIDDIT searches for discordant reads and split reads (supplementary alignments). Output vcf files are treated as intermediates and are not placed in the output folder.

#### GATK GermlineCNVCaller - CNV calling

[GATK GermlineCNVCaller](https://github.com/broadinstitute/gatk) is used to identify copy number variants in germline samples given their read counts and a model describing a sample's ploidy. Output vcf files are treated as intermediates and are not placed in the output folder.

#### CNVnator - CNV calling

[CNVnator](https://github.com/abyzovlab/CNVnator) is used to identify copy number variants in germline samples given a bam file. Output vcf files are treated as intermediates and are not placed in the output folder.

#### SVDB merge

[SVDB merge](https://github.com/J35P312/SVDB#merge) is used to merge the variant calls from GATK's GermlineCNVCaller (only if `skip_tools` doesn't include germlinecnvcaller), Manta, and TIDDIT. Output files are published in the output folder.

<details markdown="1">
<summary>Output files</summary>

- `call_sv/genome`
  - `<case_id>_sv_merge.vcf.gz`: file containing the merged variant calls. As of version 2.3.0, this file also contains mitochondrial structural variants.
  - `<case_id>_sv_merge.vcf.gz.tbi`: index of the file containing the merged variant calls.

</details>

### Variant calling - repeat expansions

#### Expansion Hunter

[Expansion Hunter](https://github.com/Illumina/ExpansionHunter) aims to estimate sizes of repeat sequences by performing a targeted search through alignments that span, flank, and are fully contained in each repeat. The files generated are ready to be used with [REViewer](https://github.com/Illumina/REViewer).

<details markdown="1">
<summary>Output files</summary>

- `repeat_expansions/`
  - `<sample_id>_repeat_expansion.vcf`: file containing variant calls.
  - `<sample_id>_exphunter_sorted.bam`: A BAMlet containing alignments of reads that overlap or located in close proximity to each variant identified by ExpansionHunter
  - `<sample_id>_exphunter_sorted.bam.bai`: Index of the BAMlet file

</details>

#### Stranger

[Stranger](https://github.com/Clinical-Genomics/stranger) annotates output files from Expansion Hunter with the pathologic implications of the repeat sizes.

<details markdown="1">
<summary>Output files</summary>

- `repeat_expansions/`
  - `<case_id>_repeat_expansion_stranger.vcf.gz`: file containing variant calls.
  - `<case_id>_repeat_expansion_stranger.vcf.gz.tbi`: index of the file containing variant calls.

</details>

### Annotation - SNV

#### bcftools roh

[bcftools roh](https://samtools.github.io/bcftools/bcftools.html#roh) is a program for detecting runs of homo/autozygosity.from only bi-allelic sites. The output files are not published in the output folder, and is passed to vcfanno for further annotation.

:::note
In the case of running a quattro, i.e. two affected children and their parents, only one of the probands will be used for annotating regions of homozygosity. This is a known limitation that we are hoping to solve in a future release.
:::

#### vcfanno

[vcfanno](https://github.com/brentp/vcfanno) allows you to quickly annotate your VCF with any number of INFO fields from any number of VCFs. It uses a simple configuration file to allow the user to specify the source annotation files and fields and how they will be added to the info of the query VCF. Values are pulled by name from the INFO field with special-cases of ID and FILTER to pull from those VCF columns. The output files are not published in the output folder, and is passed to CADD and/or VEP for further annotation.

We recommend using vcfanno to annotate SNVs with precomputed CADD scores (files can be downloaded from [here](https://cadd.gs.washington.edu/download)).

#### CADD

[CADD](https://cadd.gs.washington.edu/) is a tool for scoring the deleteriousness of single nucleotide variants as well as insertion/deletions variants in the human genome. In nf-core/raredisease, SNVs can be annotated with precomputed CADD scores using vcfanno. However, for small indels they will be calculated on the fly by CADD. The output files are not published in the output folder, and is passed to VEP for further annotation.

#### VEP

[VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) determines the effect of your variants on genes, transcripts, and protein sequence, as well as regulatory regions. We recommend annotating with the following plugins:

- LoFtool
- pLI
- SpliceAI
- MaxEntScan

Based on VEP annotations, custom scripts used by the pipeline further annotate each record with the most severe consequence, and pli scores.

> **NB**: Output files described below do not include mitochondrial annotations only if --skip_subworklows includes mt_annotation.

<details markdown="1">
<summary>Output files</summary>

- `annotate_snv/genome`
  - `<case_id>_<rhocall_vcfanno|vcfanno>_filter_<cadd_vep|vep>.vcf.gz`: file containing bcftools roh, vcfanno, cadd and vep annotations.
  - `<case_id>_<rhocall_vcfanno|vcfanno>_filter_<cadd_vep|vep>.vcf.gz.tbi`: index of the file containing bcftools roh, vcfanno, cadd and vep annotations.

</details>

#### UPD

[UPD](https://github.com/bjhall/upd) calls regions of uniparental disomy from germline exome/wgs trios. Output from UPD is passed to chromograph for making plots.

#### Chromograph

[Chromograph](https://github.com/Clinical-Genomics/chromograph) is a python package to create PNG images from genetics data such as BED and WIG files.

<details markdown="1">
<summary>Output files</summary>

- `annotate_snv/genome/*sites_chromograph`
  - `<sample_id>_rohann_vcfanno_upd_sites_<chr#>.png`: file containing a plot showing upd sites across chromosomes.
- `annotate_snv/genome/*regions_chromograph`
  - `<sample_id>_rohann_vcfanno_upd_regions_<chr#>.png`: file containing a plot showing upd regions across chromosomes.
- `annotate_snv/genome/*autozyg_chromograph`
  - `<sample_id>_rhocallviz_<chr#>.png`: file containing a plot showing regions of autozygosity across chromosomes.

</details>

#### Rhocall viz

[Rhocall viz](https://github.com/dnil/rhocall) plots binned zygosity and RHO-regions.

<details markdown="1">
<summary>Output files</summary>

- `annotate_snv/genome/<sample_id>_rhocallviz/<sample_id>_rhocallviz.bed`: file containing regions of homozygosity in bed format.
- `annotate_snv/genome/<sample_id>_rhocallviz/<sample_id>_rhocallviz.wig`: file containing the fraction of homozygous SNPs in wig format.
- `annotate_snv/genome/<sample_id>_rhocallviz/<sample_id>_rhocallviz.bw`: file containing the fraction of homozygous SNPs in bigwig format.

</details>

### Annotation - SV

#### SVDB query

[SVDB query](https://github.com/J35P312/SVDB#Query) allows you to quickly annotate your VCF with data from one or more structural variant databases. The output files are not published in the output folder, and is passed to vep for further annotation.

#### VEP

[VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) determines the effect of your variants on genes, transcripts, and protein sequence, as well as regulatory regions. We recommend annotating with pLI plugin, along with any other custom plugins you may want too use. Based on VEP annotations, custom scripts used by the pipeline further annotate each record with the most severe consequence, and pli scores.

<details markdown="1">
<summary>Output files</summary>

- `annotate_sv/`
  - `<case_id>_svdbquery_vep.vcf.gz`: file containing svdb query, and vep annotations.
  - `<case_id>_svdbquery_vep.vcf.gz.tbi`: index of the file containing bcftools roh, vcfanno, and vep annotations.

</details>

### Mitochondrial analysis

Mitochondrial analysis is run by default. If you want to turn off annotations set `--skip_subworkflows mt_annotation`.

#### Haplocheck

[Haplocheck](https://github.com/genepi/haplocheck/) analyses the mitochondrial content to detect contamination in samples. The results are displayed in MultiQC.

<details markdown="1">
<summary>Output files</summary>

- `{outputdir}/haplocheck/`
  - `<sampleid>.raw.txt`: haplocheck raw output. Read more about the file content [here](https://mitoverse.readthedocs.io/haplocheck/haplocheck/#textual-report-file-raw).

</details>

#### Alignment and variant calling

[Alignment and variant calling - GATK Mitochondrial short variant discovery pipeline](https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels-) The mitochondrial genome poses several challenges to the identification and understanding of somatic variants. The circularity of the mitochondrial genome means that the breakpoint in the reference genome is at an arbitrary position in the non-coding control region, creating a challenge in analyzing variation. Additionally, insertions of mitochondrial DNA into the nuclear genome (NuMTs) complicate the mapping of the mitochondrial genome and the distinction between NuMTs and the mitochondrial contig of interest. Lastly, mitochondrial variants often have very low heteroplasmy. Such low allele fraction (AF) variants can thus be mistaken for inherent sequencer noise.

The pipeline for mitochondrial variant discovery, using Mutect2, uses a high sensitivity to low AF and separate alignments using opposite genome breakpoints to allow for the tracing of lineages of rare mitochondrial variants.

- `call_snv/mitochondria`
  - `<case_id>_mitochondria.vcf.gz`: normalized vcf file containing MT variants.
  - `<case_id>_mitochondria.vcf.gz.tbi`: index of the vcf file containing MT variants.

##### MT deletion script

[MT deletion script](https://github.com/dnil/mitosign/blob/master/run_mt_del_check.sh) lists the fraction of mitochondrially aligning read pairs (per 1000) that appear discordant, as defined by an insert size of more than 1.2 kb (and less than 15 kb due to the circular nature of the genome) using samtools.

- `call_sv/mitochondria`
  - `<sample_id>_mitochondria_deletions.txt`: file containing deletions.

##### eKLIPse

[eKLIPse](https://github.com/dooguypapua/eKLIPse) allows the detection and quantification of large mtDNA rearrangements.

- `call_sv/mitochondria`
  - `eKLIPse_<sample_id>_deletions.csv`: file containing all predicted deletions.
  - `eKLIPse_<sample_id>_genes.csv`: file summarizing cumulated deletions per mtDNA gene.
  - `eKLIPse_<sample_id>.png`: circos plot.

#### Annotation

##### HaploGrep2

[HaploGrep2](https://github.com/seppinho/haplogrep-cmd) allows detecting artificial recombinants and missing variants as well as annotating rare and phantom mutations in mitochondria. Haplogrep generates a text report, which is published by default.

<details markdown="1">
<summary>Output files</summary>

- `annotate_snv/mitochondria`
  - `<case_id>*haplogrep.txt`: file containing haplogroup information.

</details>

##### vcfanno

[vcfanno](https://github.com/brentp/vcfanno) allows you to quickly annotate your VCF with any number of INFO fields from any number of VCFs. It uses a simple conf file to allow the user to specify the source annotation files and fields and how they will be added to the info of the query VCF. Values are pulled by name from the INFO field with special-cases of ID and FILTER to pull from those VCF columns. The output files are not published in the output folder, and is passed to vep for further annotation.

We recommend using vcfanno to annotate SNVs with precomputed CADD scores (files can be downloaded from [here](https://cadd.gs.washington.edu/download)).

##### CADD

[CADD](https://cadd.gs.washington.edu/) is a tool for scoring the deleteriousness of single nucleotide variants as well as insertion/deletions variants in the human genome. In nf-core/raredisease, SNVs can be annotated with precomputed CADD scores using vcfanno. However, for small indels they will be calculated on the fly by CADD. The output files are not published in the output folder, and is passed to VEP for further annotation.

##### Hmtnote

[HmtNote](https://github.com/robertopreste/HmtNote) annotates vcf containing human mitochondrial variants with HmtVar. It will run offline by default with a database within the container.

##### VEP

[VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) determines the effect of your variants on genes, transcripts, and protein sequence, as well as regulatory regions.

<details markdown="1">
<summary>Output files</summary>

- `annotate_snv/mitochondria`
  - `<case_id>_mitochondria_hmtnote_vcfanno_<cadd_vep|vep>.vcf.gz`: file containing mitochondrial annotations.
  - `<case_id>_mitochondria_hmtnote_vcfanno_<cadd_vep|vep>.vcf.gz.tbi`: index of the file containing mitochondrial annotations.

</details>

### Filtering and ranking

#### Filter_vep

[filter_vep from VEP](https://www.ensembl.org/info/docs/tools/vep/script/vep_filter.html) is used to subset the variants based on a list of HGNC ID:s. Typical use case is that you want to filter your results to only include variants in a predefined set of clinically relevant genes. This step is optional and can be disabled by using the flag `--skip_subworkflows generate_clinical_set`. You will always get the complete VCF together with the clinical VCF.

#### GENMOD

[GENMOD](https://github.com/Clinical-Genomics/genmod) is a simple to use command line tool for annotating and analyzing genomic variations in the VCF file format. GENMOD can annotate genetic patterns of inheritance in vcf files with single or multiple families of arbitrary size. Each variant will be assigned a predicted pathogenicity score. The score will be given both as a raw score and a normalized score with values between 0 and 1. The tags in the INFO field are `RankScore` and `RankScoreNormalized`. The score can be configured to fit your annotations and preferences by modifying the score config file.

<details markdown="1">
<summary>Output files</summary>

- `rank_and_filter/`
  - `<case_id>_mt_ranked_clinical.vcf.gz`: file containing clinically relevant mitochondrial SNVs, and only contains variants less than 5%VAF by default.
  - `<case_id>_mt_ranked_clinical.vcf.gz.tbi`: index of the file containing clinically relevant mitochondrial SNVs.
  - `<case_id>_mt_ranked_research.vcf.gz`: file containing mitochondrial SNV annotations with their rank scores.
  - `<case_id>_mt_ranked_research.vcf.gz.tbi`: index of the file containing mitochondrial SNV annotations with their rank scores.
  - `<case_id>_snv_ranked_clinical.vcf.gz`: file containing clinically relevant SNVs (does not include mitochondrial variants).
  - `<case_id>_snv_ranked_clinical.vcf.gz.tbi`: index of the file containing clinically relevant SNVs.
  - `<case_id>_snv_ranked_research.vcf.gz`: file containing SNV annotations with their rank scores (does not include mitochondrial variants).
  - `<case_id>_snv_ranked_research.vcf.gz.tbi`: index of the file containing SNV annotations with their rank scores.
  - `<case_id>_sv_ranked_clinical.vcf.gz`: file containing clinically relevant SVs (includes mitochondrial variants).
  - `<case_id>_sv_ranked_clinical.vcf.gz.tbi`: index of the file containing clinically relevant SVs.
  - `<case_id>_sv_ranked_research.vcf.gz`: file containing SV annotations with their rank scores (includes mitochondrial variants).
  - `<case_id>_sv_ranked_resarch.vcf.gz.tbi`: index of the file containing SV annotations with their rank scores.

</details>

### Mobile element analysis

#### Calling mobile elements

Mobile elements are identified from the bam file using [RetroSeq](https://github.com/tk2/RetroSeq) and the indiviual calls are merged to case VCF using SVDB.

<details markdown="1">
<summary>Output files</summary>

- `call_mobile_elements/`
  - `<case_id>_mobile_elements.vcf.gz`: file containing mobile elements.
  - `<case_id>_mobile_elements.vcf.gz.tbi`: index of the file containing mobile elements.

</details>

#### Annotating mobile elements

The mobile elements are annotated with allele frequencies and allele counts using SVDB. These annotation files needed are preferably produced from a representative population. Further annoation is done using VEP and the resulting VCF is filtered using bcftools. The default filter is to only keep elements with `PASS` in the filter column but if no other post-processing is done we reccomend supplementing with an exclude expression based on population allele frequencies. The filtering key is dependent on the annotation files used but an example expression could look like this: `--exclude 'INFO/swegen_sva_FRQ > 0.1'`. If a list of HGNC id's have been supplied with the option `--vep_filters`, variants matching those id's will be presented in a seperate file using [filter_vep from VEP](https://www.ensembl.org/info/docs/tools/vep/script/vep_filter.html). This option can be disabled using the flag `--skip_subworkflows generate_clinical_set`. A VCF corresponding to the complete set of variants will also be produced.

<details markdown="1">
<summary>Output files</summary>

- `rank_and_filter/`
  - `<case_id>_mobile_elements_annotated_research.vcf.gz`: VCF containting the complete set of annotated mobile elements.
  - `<case_id>_mobile_elements_annotated_research.vcf.gz.tbi`: Index for VCF containting the complete set of annotated mobile elements.
  - `<case_id>_mobile_elements_annotated_clinical.vcf.gz`: VCF containing selected annotated mobile elements.
  - `<case_id>_mobile_elements_annotated_clincial.vcf.gz.tbi`: Index for VCF containing selected annotated mobile elements.

</details>

### Variant evaluation

Provided a truth set, SNVs can be evaluated using RTG Tools' vcfeval engine. Output files generated are listed below with a short description, but if you'd like to know more about what's in each of the files, refer to RTG Tools documentation [here](https://www.animalgenome.org/bioinfo/resources/manuals/RTGOperationsManual.pdf).

<details markdown="1">
<summary>Output files</summary>

- `rtgvcfeval/`

  - `<sample_id>_vcfeval.fn.vcf.gz`: contains variants from the baseline VCF which were not correctly called.
  - `<sample_id>_vcfeval.fn.vcf.gz.tbi`: index of the \*fn.vcf file
  - `<sample_id>_vcfeval.fp.vcf.gz`: contains variants from the calls VCF which do not agree with baseline variants.
  - `<sample_id>_vcfeval.fp.vcf.gz.tbi`: index of the \*fp.vcf file
  - `<sample_id>_vcfeval.non_snp_roc.tsv.gz`: contains ROC data derived from those variants which were not represented as
    SNPs.
  - `<sample_id>_vcfeval.phasing.txt`: containing the data on the phasing
  - `<sample_id>_vcfeval.snp_roc.tsv.gz`: contains ROC data derived from only those variants which were represented as SNPs.
  - `<sample_id>_vcfeval.summary.txt`: contains the match summary statistics printed to standard output.
  - `<sample_id>_vcfeval.tp-baseline.vcf.gz`: contains those variants from the baseline VCF which agree with variants in the
    calls VCF.
  - `<sample_id>_vcfeval.tp-baseline.vcf.gz.tbi`: index of the \*tp-baseline.vcf file
  - `<sample_id>_vcfeval.tp.vcf.gz`: contains those variants from the calls VCF which agree with variants in the baseline VCF
  - `<sample_id>_vcfeval.tp.vcf.gz.tbi`: index of the \*tp.vcf file
  - `<sample_id>_vcfeval.weighted_roc.tsv.gz`: contains ROC data derived from all analyzed call variants, regardless of their
    representation.

</details>

### Gens

The sequencing data can be prepared for visualization of CNVs in [Gens](https://github.com/Clinical-Genomics-Lund/gens). You can turn it off by supplying the option `--skip_tools gens`. You can read more about how to setup Gens [here](https://github.com/Clinical-Genomics-Lund/gens).

<details markdown="1">
<summary>Output files</summary>

- `gens/`

  - `<sample_id>_gens.baf.bed.gz`: contains sample b-allele frequencies in bed format.
  - `<sample_id>_gens.baf.bed.gz.tbi`: index of the \*baf.bed.gz file.
  - `<sample_id>_gens.cov.bed.gz`: contains sample coverage in bed format.
  - `<sample_id>_gens.cov.bed.gz.tbi`: index of the \*cov.bed.gz file.

</details>

### Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>
