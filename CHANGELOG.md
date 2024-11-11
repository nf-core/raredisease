# nf-core/raredisease: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 2.3.0dev - Getafix [xxxx-xx-xx]

### `Added`

- A new analysis option `mito` to call and annotate only mitochondrial variants [#608](https://github.com/nf-core/raredisease/pull/608)

### `Changed`

- Suffix used to identify unique fastq pairs from "\_T" to "\_LNUMBER" [#638](https://github.com/nf-core/raredisease/pull/638)
- Merge output from germlinecnvcaller [#635](https://github.com/nf-core/raredisease/pull/635)
- Update tools [#623](https://github.com/nf-core/raredisease/pull/623)
- Update output file name prefix for upd and chromograph to sample-based [#620](https://github.com/nf-core/raredisease/pull/620)
- Update tools [#619](https://github.com/nf-core/raredisease/pull/619)
- Report only variants above 5% heteroplasmy in the clinical vcf file for mitochondria [#616](https://github.com/nf-core/raredisease/pull/616)

### `Fixed`

- Restrict deepvariant analysis of WES samples to bait regions [#633](https://github.com/nf-core/raredisease/pull/633)
- bcftools annotate declaration in annotate CADD subworkflow [#624](https://github.com/nf-core/raredisease/pull/624)
- Rhocallviz subworkflow will only be invocated once per sample [#621](https://github.com/nf-core/raredisease/pull/621)
- Allow for VEP version 112 to be used and set it to default [#617](https://github.com/nf-core/raredisease/pull/617)
- Updated createCaseChannel function to include a check for maternal and paternal ids being set to a numeric 0 [#643](https://github.com/nf-core/raredisease/pull/643)

### Parameters

### Tool updates

| Tool       | Old version | New version |
| ---------- | ----------- | ----------- |
| bcftools   | 1.18        | 1.20        |
| ensemblvep | 112         | 113         |
| genmod     | 3.8.2       | 3.9         |
| mosdepth   | 0.3.6       | 0.3.8       |
| multiqc    | 1.21        | 1.25.1      |
| picard     | 3.1.1       | 3.3.0       |
| samtools   | 1.19.2      | 1.21        |
| sentieon   | 202308.02   | 202308.03   |
| stranger   | 0.8.1       | 0.9.2       |
| svdb       | 2.8.1       | 2.8.2       |
| tabix      | 1.19.1      | 1.20        |

## 2.2.0 - Dogmatix [2024-09-13]

### `Added`

- A new parameter `mt_aligner` to control which aligner is used to align reads to mitochondria [#600](https://github.com/nf-core/raredisease/pull/600)
- A new parameter `par_bed` to pass a PAR bed files to deepvariant [#598](https://github.com/nf-core/raredisease/pull/598)
- A new functionality to pass gzipped resources to vcfanno_extra_resources [#589](https://github.com/nf-core/raredisease/pull/589)
- A new parameter `vcfanno_extra_resources` to pass an extra resource to vcfanno [#588](https://github.com/nf-core/raredisease/pull/588)
- A new parameter `scatter_count` to control how many interval files are created from a genome (used to parallelize annotations) [#585](https://github.com/nf-core/raredisease/pull/585)
- Print warning messages if user intends to perform ranking when there are no affected samples [#579](https://github.com/nf-core/raredisease/pull/579)
- Two new parameters `skip_repeat_annotation` and `skip_repeat_calling` to skip calling and annotation of repeat expansions [#574](https://github.com/nf-core/raredisease/pull/574)
- A new parameter `skip_smncopynumbercaller` to skip smncopynumbercaller module [#574](https://github.com/nf-core/raredisease/pull/574)
- A new parameter `skip_sv_calling` to skip sv calling workflow [#572](https://github.com/nf-core/raredisease/pull/572)
- Two new parameters `skip_snv_calling` and `skip_repeat_analysis` to skip snv calling and repeat analysis respectively [#571](https://github.com/nf-core/raredisease/pull/571)
- Two new parameters `mbuffer_mem` and `samtools_sort_threads` to control resources given to mbuffer and samtools sort in the bwameme module [#570](https://github.com/nf-core/raredisease/pull/570)

### `Changed`

- Update default vep container from v110-v112 [#609](https://github.com/nf-core/raredisease/pull/609)
- Default index for vcfanno extra annotation files from tbi to csi [#606](https://github.com/nf-core/raredisease/pull/606)
- Updated the model for Sentieon DNAScope to v1.1 [#601](https://github.com/nf-core/raredisease/pull/601)
- bwameme can no longer be used to align mitochondrial reads [#600](https://github.com/nf-core/raredisease/pull/600)
- Males' X and Y chromosomes will be treated as haploids during variant calling by deepvariant [#598](https://github.com/nf-core/raredisease/pull/598)
- Acceptable type for lane field in the samplesheet from number to string [#597](https://github.com/nf-core/raredisease/pull/597)
- Allow `0` as a valid value for `sex` in the samplesheet [#595](https://github.com/nf-core/raredisease/pull/595)
- Updated deepvariant to version 1.6.1 [#587](https://github.com/nf-core/raredisease/pull/587)
- Parallelized vcfanno [#585](https://github.com/nf-core/raredisease/pull/585)
- Skip ROH calling with bcftools if there are no affected samples [#579](https://github.com/nf-core/raredisease/pull/579)
- Refactored tool citation list [#577](https://github.com/nf-core/raredisease/pull/577)
- Removed `skip_repeat_analysis` added in #571 [#574](https://github.com/nf-core/raredisease/pull/574)
- Remove several skip parameters that had been included in the pipeline to avoid failed CI tests (see parameters table below) [#574](https://github.com/nf-core/raredisease/pull/574)
- `readcount_intervals` parameter is now mandatory for running germlinecnvcaller. [#570](https://github.com/nf-core/raredisease/pull/570)
- Turn off CNVnator, TIDDIT, SMNCopyNumberCaller, Gens, and Vcf2cytosure for targeted analysis [#573](https://github.com/nf-core/raredisease/pull/573)

### `Fixed`

- Issues that cropped up when `aligner` and `mt_aligner` were different [#605](https://github.com/nf-core/raredisease/pull/605)
- Update docs to show 'vep_plugin_files' as a mandatory parameter for SNV annotation [#594](https://github.com/nf-core/raredisease/pull/594)
- Error in SVDB merge when only a single SV caller is run [#586](https://github.com/nf-core/raredisease/pull/586)
- Errors due to misplaced version statements [#578](https://github.com/nf-core/raredisease/pull/578)
- Stub crashes due to peddy reported in [#566](https://github.com/nf-core/raredisease/issues/566) [#576](https://github.com/nf-core/raredisease/pull/576]
- Docker manifest error from gnu-wget container [#570](https://github.com/nf-core/raredisease/pull/570)
- Citations for bwameme [#563](https://github.com/nf-core/raredisease/pull/563)

### Parameters

| Old parameter   | New parameter            |
| --------------- | ------------------------ |
|                 | mbuffer_mem              |
|                 | mt_aligner               |
|                 | samtools_sort_threads    |
|                 | skip_repeat_calling      |
|                 | skip_snv_calling         |
|                 | skip_sv_calling          |
| skip_eklipse    |                          |
| skip_fastqc     |                          |
| skip_haplocheck |                          |
| skip_qualimap   |                          |
|                 | skip_smncopynumbercaller |
|                 | skip_repeat_annotation   |
|                 | scatter_count            |
|                 | vcfanno_extra_resources  |

### Tool updates

| Tool        | Old version | New version |
| ----------- | ----------- | ----------- |
| Deepvariant | 1.5.0       | 1.6.1       |
| ensemblvep  | 110         | 112         |

## 2.1.0 - Obelix [2024-05-29]

### `Added`

- A new aligner, bwameme [#553](https://github.com/nf-core/raredisease/pull/553)
- A new parameter `run_mt_for_wes` to turn on mitochondrial analysis for targeted analysis [#552](https://github.com/nf-core/raredisease/pull/552)
- A new parameter `bwa_as_fallback` to switch aligner to bwa in case bwamem2 fails [#551](https://github.com/nf-core/raredisease/pull/551)
- A new parameter `skip_me_calling` to skip mobile element calling and the subsequent annotation of them [#556](https://github.com/nf-core/raredisease/pull/556)

### `Changed`

- Changed valid values for sex according to the PED file format [#550](https://github.com/nf-core/raredisease/pull/550)
- Refactored config files [#538](https://github.com/nf-core/raredisease/pull/538)
- Refactored mobile element annotation subworkflow files [#538](https://github.com/nf-core/raredisease/pull/538)
- Refactored to remove "a process is defined more than once" warning [#557](https://github.com/nf-core/raredisease/pull/557)
- Updated modules [#558](https://github.com/nf-core/raredisease/pull/558)

### `Fixed`

- Include multiallelic indel sites in CADD scoring jobs [#545](https://github.com/nf-core/raredisease/pull/545)
- Fixed issues with samtools merge not being run on samples sequenced over multiple lanes [#538](https://github.com/nf-core/raredisease/pull/538)
- Fixed join issues in the mobile element calling subworkflow which occured when mobile_element_references were not provided [#556](https://github.com/nf-core/raredisease/pull/556)

### Parameters

| Old parameter | New parameter   |
| ------------- | --------------- |
|               | bwameme         |
|               | bwa_as_fallback |
|               | run_mt_for_wes  |
|               | skip_me_calling |

:::note
Parameter has been updated if both old and new parameter information is present.
Parameter has been added if just the new parameter information is present.
Parameter has been removed if new parameter information isn't present.
:::

### Module updates

| Tool     | Old version | New version |
| -------- | ----------- | ----------- |
| bwa      | 0.7.17      | 0.7.18      |
| CADD     | 1.6.1       | 1.6.post1   |
| Sentieon | 202308.01   | 202308.02   |
| bwameme  |             | 1.0.6       |

:::note
Version has been updated if both old and new version information is present.
Version has been added if just the new version information is present.
Version has been removed if new version information isn't present.
:::

## 2.0.1 - Asterix (Patch) [2024-03-25]

### `Fixed`

- Germlinecnvcaller subworkflow uses the output channel `casecalls` from germlinecnvcaller module instead of `calls` which was invalid. [#535](https://github.com/nf-core/raredisease/issues/535)

## 2.0.0 - Asterix [2024-03-18]

### `Added`

- Use `nf-validation` plugin for parameter and samplesheet validation [#386](https://github.com/nf-core/raredisease/pull/386)
- A new parameter `skip_vep_filter` to skip filtering based on vep results [#416](https://github.com/nf-core/raredisease/pull/416)
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
- Added a new parameter `vep_filters_scout_fmt` to supply a bed-like file exported by scout to be used in filter_vep [#511](https://github.com/nf-core/raredisease/pull/511).
- Added two new parameters `variant_consequences_snv` and `variant_consequences_sv` to supply variant consequence files for annotating SNVs and SVs. [#509](https://github.com/nf-core/raredisease/pull/509)

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
- Replaced the logic where we added an arbitrary substring to keep file names unique after alignment which we then removed using a split operator, with a simple copy operation. [#425](https://github.com/nf-core/raredisease/pull/425)
- Preventing a crash of rhocall annotate in the case of running four individuals whereof two are affected.
- Fixed memory qualifier in gatk4 germlinecnvcaller and postprocessgermlinecnvcalls
- Fixed wrong process names when outputting versions in `ALIGN_SENTIEON` and `CALL_SNV`.
- Fixed gens subworkflow [#515](https://github.com/nf-core/raredisease/pull/515)

### Parameters

| Old parameter         | New parameter                         |
| --------------------- | ------------------------------------- |
|                       | `--cnvnator_binsize`                  |
|                       | `--gens_pon_female`                   |
|                       | `--gens_pon_male`                     |
|                       | `--min_trimmed_length`                |
|                       | `--mobile_element_references`         |
|                       | `--mobile_element_svdb_annotations`   |
|                       | `--mt_subsample_rd`                   |
|                       | `--mt_subsample_seed`                 |
|                       | `--ngsbits_samplegender_method`       |
|                       | `--rtg_truthvcfs`                     |
|                       | `--run_rtgvcfeval`                    |
|                       | `--sample_id_map`                     |
|                       | `--score_config_mt`                   |
|                       | `--sdf`                               |
| `--pcr_amplification` | `--sentieon_dnascope_pcr_indel_model` |
|                       | `--skip_eklipse`                      |
|                       | `--skip_fastqc`                       |
|                       | `--skip_fastp`                        |
|                       | `--skip_gens`                         |
|                       | `--skip_germlinecnvcaller`            |
|                       | `--skip_haplocheck`                   |
|                       | `--skip_me_annotation`                |
|                       | `--skip_mt_annotation`                |
|                       | `--skip_mt_subsample`                 |
|                       | `--skip_peddy`                        |
|                       | `--skip_qualimap`                     |
|                       | `--skip_vcf2cytosure`                 |
|                       | `--skip_vep_filter`                   |
|                       | `--svdb_query_bedpedbs`               |
|                       | `--variant_consequences_snv`          |
|                       | `--variant_consequences_sv`           |
|                       | `--vcf2cytosure_blacklist`            |
|                       | `--vep_plugin_files`                  |
|                       | `--vep_filters_scout_fmt`             |
| `--gens_pon`          |                                       |
| `--gens_switch`       |                                       |
| `--skip_cnv_calling`  |                                       |
| `--skip_mt_analysis`  |                                       |

:::note
Parameter has been updated if both old and new parameter information is present.
Parameter has been added if just the new parameter information is present.
Parameter has been removed if new parameter information isn't present.
:::

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
