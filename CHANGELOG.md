# nf-core/raredisease: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0 - Genie [xxxx-xx-xx]

### `Added`

- Use `nf-validation` plugin for parameter and samplesheet validation [#386](https://github.com/nf-core/raredisease/pull/386)
- A new parameter to skip filtering based on vep results [#416](https://github.com/nf-core/raredisease/pull/416)
- A `metromap` representating the core parts of the pipeline [#428](https://github.com/nf-core/raredisease/pull/428)
- Metromap and logos for light and dark theme [#432](https://github.com/nf-core/raredisease/pull/432)
- New parameters to skip qualimap and eklipse (`--skip_qualimap` and `--skip_eklipse`) [#436](https://github.com/nf-core/raredisease/pull/436)
- Fix "there is no process matching config selector warnings" [#435](https://github.com/nf-core/raredisease/pull/435)
- New parameters to skip fastqc and haplocheck (`--skip_fastqc` and `--skip_haplocheck`) [#438](https://github.com/nf-core/raredisease/pull/438)
- CNVnator for copy number variant calling [#438](https://github.com/nf-core/raredisease/pull/434)
- A new parameter `svdb_query_bedpedbs` to provide bedpe files as databases for SVDB query [#449](https://github.com/nf-core/raredisease/pull/449)
- ngsbits samplegender to check sex [#453](https://github.com/nf-core/raredisease/pull/453)
- New workflow for generating cgh files from SV vcfs for interpretation in the CytosSure interpretation software. Turned off by default [#456](https://github.com/nf-core/raredisease/pull/456/)
- Fastp to do adapter trimming. It can be skipped using `--skip_fastp` [#457](https://github.com/nf-core/raredisease/pull/457)
- New workflow for calling insertion of mobile elements [#440](https://github.com/nf-core/raredisease/pull/440)
- GATK CNVCaller uses segments instead of intervals, filters out "reference" segments between the calls, and fixes a bug with how `ch_readcount_intervals` was handled [#472](https://github.com/nf-core/raredisease/pull/472)
- bwa aligner [#474](https://github.com/nf-core/raredisease/pull/474)
- Add FOUND_IN tag, which mentions the variant caller that found the mutation, in the INFO column of the vcf files [#471](https://github.com/nf-core/raredisease/pull/471)
- A new parameter `vep_plugin_files` to supply files required by vep plugins [#482](https://github.com/nf-core/raredisease/pull/482)
- New workflow for annotating mobile elements [#483](https://github.com/nf-core/raredisease/pull/483)
- Added a functionality to subsample mitochondrial alignment, and a new parameter `skip_mt_subsample` to skip the subworkflow [#508](https://github.com/nf-core/raredisease/pull/508).
- Chromograph to plot coverage across chromosomes [#507](https://github.com/nf-core/raredisease/pull/507)

### `Changed`

- Tiddit updated to v3.6.1 [#385](https://github.com/nf-core/raredisease/pull/385)
- Installed the nf-core version of the sentieon/bwamemindex module [#397](https://github.com/nf-core/raredisease/pull/397)
- Installed the nf-core version of the sentieon/bwamem module [#398](https://github.com/nf-core/raredisease/pull/398)
- Installed the nf-core version of the sentieon/readwriter module [#399](https://github.com/nf-core/raredisease/pull/399)
- Installed the nf-core version of the sentieon/datametrics module [#400](https://github.com/nf-core/raredisease/pull/400)
- Installed the nf-core version of the sentieon/dedup module. The dedup module also contains a call to Sentieon's LocusCollector [#401](https://github.com/nf-core/raredisease/pull/401)
- Removing Sentieon-based BQSR. Recent Illumina sequencers tend to provide well-calibrated BQs, so BQSR may not provide much benefit [#402](https://github.com/nf-core/raredisease/pull/402)
- Installed the nf-core version of the sentieon/dnamodelapply module [#403](https://github.com/nf-core/raredisease/pull/403)
- Installed the nf-core version of the sentieon/wgsmetricsalgo module [#404](https://github.com/nf-core/raredisease/pull/404)
- Installed the nf-core version of the sentieon/dnascope module [#406](https://github.com/nf-core/raredisease/pull/406)
- Breaks down mitochondrial analysis workflow into smaller subworkflows that are more modular [#419](https://github.com/nf-core/raredisease/pull/419)
- Replaced the parameter skip_mt_analysis which was used to turn on/off the mitochondrial workflow [#419](https://github.com/nf-core/raredisease/pull/419)
- Adds a new parameter skip_mt_annotation which can be used to turn on/off annotation and ranking for mitochondrial SNVs [#419](https://github.com/nf-core/raredisease/pull/419)
- Changed the name of the parameter from `skip_cnv_calling` to `skip_germlinecnvcaller` [#435](https://github.com/nf-core/raredisease/pull/435)
- Check SVDB query input files for existence and correct format [#476](https://github.com/nf-core/raredisease/pull/476)
- Change hardcoded platform value to params.platform in align_MT.config [#475](https://github.com/nf-core/raredisease/pull/475)
- The split into clincial and research VCFs is now done before ranking the varaints [#485](https://github.com/nf-core/raredisease/pull/485)
- Installed the nf-core version of ensemblvep/vep module [#482](https://github.com/nf-core/raredisease/pull/482)
- The filenames of the ranked output VCF files have been changed. See [output.md](docs/output.md#filtering-and-ranking) for more information[#485](https://github.com/nf-core/raredisease/pull/485)
- Patched cnvnator module so that the processes didn't have to rerun after a failed run [#503](https://github.com/nf-core/raredisease/pull/503).
- Added a local module to generate bed files with variant caller names [#505](https://github.com/nf-core/raredisease/pull/505).

### `Fixed`

- Invalid GATK4 container which caused incorrect singularity downloads with nf-core download [nf-core/modules #3668](https://github.com/nf-core/modules/issues/3668)
- Make the default cram prefix same as markduplicates prefix [#392](https://github.com/nf-core/raredisease/pull/392)
- Sort ranked SV vcf before indexing with tabix [#393](https://github.com/nf-core/raredisease/pull/393)
- Make target bed file optional for WGS mode (Issue [#375](https://github.com/nf-core/raredisease/issues/375)) [#395](https://github.com/nf-core/raredisease/pull/395)
- Added constraints to block the pipeline from running CollectWgsMetrics on WES samples [#396](https://github.com/nf-core/raredisease/pull/396)
- Updated modules from nf-core [#412](https://github.com/nf-core/raredisease/pull/412)
- If present, remove duplicate entries in probands and upd_children in the meta. [#420](https://github.com/nf-core/raredisease/pull/420)
- Fixes vep starting as many instances as the square of the number of scatters. [#405](https://github.com/nf-core/raredisease/pull/405)
- Replaced the logic where we added an arbitrary substring to keep file names unique after alignment which we then removed using a split operator, with a simple copy operation. [#425](https://github.com/nf-core/raredisease/pull/425/files)
- Preventing a crash of rhocall annotate in the case of running four individuals whereof two are affected.
- Fixed memory qualifier in gatk4 germlinecnvcaller and postprocessgermlinecnvcalls
- Fixed wrong process names when outputting versions in `ALIGN_SENTIEON` and `CALL_SNV`.

### `Updated`

| Tool         | Old version | New version |
| ------------ | ----------- | ----------- |
| `filter_vep` | 107         | 110         |
| `multiqc`    | 1.14        | 1.15        |
| `untar`      | (grep 3.4)  | (grep 3.11) |

## v1.1.1 - Abu (Patch) [2023-07-26]

### `Fixed`

- Avoids errors thrown by bcftools concat due to sample names in input vcf files not being in same order [#388](https://github.com/nf-core/raredisease/pull/388)

## v1.1.0 - Abu [2023-07-21]

### `Added`

- Add GATK's cnv calling pipeline [#362](https://github.com/nf-core/raredisease/pull/362)
- GATK's ShiftFasta to generate all the files required for mitochondrial analysis [#354](https://github.com/nf-core/raredisease/pull/354)
- Feature to calculate CADD scores for indels [#325](https://github.com/nf-core/raredisease/pull/325)
- HmtNote to annotate mitochondria [#355](https://github.com/nf-core/raredisease/pull/355)
- MT del script to detect mitochondrial deletions [#349](https://github.com/nf-core/raredisease/pull/349)
- eKLIPse to identify large mitochondrial deletions [#365](https://github.com/nf-core/raredisease/pull/365)
- UPD+Chromograph to identify and visualize UPD sites and regions in the chromosomes [#364](https://github.com/nf-core/raredisease/pull/364) and [#366](https://github.com/nf-core/raredisease/pull/366)
- Added check for presence of case id for each sample in samplesheet [#357](https://github.com/nf-core/raredisease/pull/357)

### Fixed

- Avoiding publishing uncompressed VCF-file from `HMTNOTE_ANNOTATE`. (The corresponding compressed VCF-file still gets published.) [#368](https://github.com/nf-core/raredisease/pull/368)

## v1.0.0 - Aladdin [2023-06-01]

Initial release of nf-core/raredisease, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- FastQC read quality control
- Read mapping with BWAmem2/Sentieon
- Qualimap & Picard tools quality control metrics
- Call repeat expansions with ExpansionHunter and Stranger
- SNV calling with DeepVariant/Sentieon
- SV calling with Manta and TIDDIT
- SNV annotation with bcftools roh, vcfanno, and vep
- SV annotation with SVDB query and vep
- Separate workflow for analysing and annotating mitochondrial variants
- Call copy number variants in the SMN gene using SMNCopyNumberCaller
