# nf-core/raredisease: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 3.2.0dev - Luigi [XXXX-XX-XX]

### `Added`

- Add a real (non-stub) test to `annotate_genome_snvs` using the new minimal 9-region GIAB dataset in `--gtf` VEP mode, and migrate `annotate_mt_snvs`'s existing real test from the old large-genome fixtures to the same minimal dataset [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #984](https://github.com/nf-core/raredisease/pull/984)
- Add a real (non-stub) test to `generate_cytosure_files` using the new minimal 9-region GIAB dataset, and migrate `call_sv_manta`'s existing real test to the same minimal dataset [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #986](https://github.com/nf-core/raredisease/pull/986)
- Add a real (non-stub) test to `annotate_rhocallviz` using the new minimal 9-region GIAB dataset, and migrate `annotate_consequence_pli`'s real tests off the generic nf-core/modules VCF fixture onto a genuine VEP-CSQ-annotated VCF derived from the same minimal dataset [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #985](https://github.com/nf-core/raredisease/pull/985)
- Add a real (non-stub) test to `annotate_structural_variants` using the new minimal 9-region GIAB dataset, and migrate `call_sv`'s existing real test from the old large-genome fixtures to the same minimal dataset [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #983](https://github.com/nf-core/raredisease/pull/983)
- Add `somalier/extract` and `somalier/relate` modules via the `nf-core/vcf_extract_relate_somalier` subworkflow, ported from `patch` [issue #446](https://github.com/nf-core/raredisease/issues/446) [PR #891](https://github.com/nf-core/raredisease/pull/891)
- Add the `fastdup` module as a faster alternative to Picard MarkDuplicates in `align_bwa_bwamem2_bwameme`, ported from `patch` [issue #864](https://github.com/nf-core/raredisease/issues/864) [PR #876](https://github.com/nf-core/raredisease/pull/876)
- Add a real (non-stub) test to `call_snv_sentieon`, covering the two-sample merge scenario (`BCFTOOLS_MERGE`) to match its stub test's coverage, and migrate `call_snv_deepvariant`'s real test to the new minimal 9-region GIAB dataset [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #987](https://github.com/nf-core/raredisease/pull/987)
- Add `test_vcf_singleton` profile [issue #869](https://github.com/nf-core/raredisease/issues/869) [PR #978](https://github.com/nf-core/raredisease/pull/978)
- Add `test_align_singleton` profile (BAM-only singleton input) [issue #869](https://github.com/nf-core/raredisease/issues/869) [PR #977](https://github.com/nf-core/raredisease/pull/977)
- Add `peddy_sites` parameter: lets a custom sites file (`chrom:pos:ref:alt` per line) be used in place of peddy's bundled hg19/hg38 site lists, which are absolute-coordinate and won't overlap a coordinate-remapped or sliced reference, causing `het_check` to crash on zero overlapping sites [issue #869](https://github.com/nf-core/raredisease/issues/869) [PR #971](https://github.com/nf-core/raredisease/pull/971)
- Add GTF-mode VEP annotation support: `ENSEMBLVEP_VEP` can now annotate directly from a bgzipped, tabix-indexed GTF file (`vep_gtf`/`vep_gtf_tbi`) instead of an Ensembl cache, for custom/non-standard references with no corresponding Ensembl cache (e.g. a sliced or coordinate-remapped genome); mutually exclusive with `vep_cache` [issue #869](https://github.com/nf-core/raredisease/issues/869) [PR #970](https://github.com/nf-core/raredisease/pull/970)
- Added `pre_vep_snv_filter_expression` parameter to configure the `bcftools view --exclude` expression used to filter SNVs before VEP annotation, instead of hardcoding it in `conf/modules/annotate_genome_snvs.config` [issue #929](https://github.com/nf-core/raredisease/issues/929) [PR #958](https://github.com/nf-core/raredisease/pull/958)
- Add a VCF entry point: samplesheets can now supply precalled, case-level SNV/SV/MT VCFs directly, skipping the corresponding calling step and feeding straight into annotation and ranking. See `docs/usage.md` for samplesheet details [issue #261](https://github.com/nf-core/raredisease/issues/261) [PR #936](https://github.com/nf-core/raredisease/pull/936)
- Add test coverage for the VCF entry point: function-level nf-test cases for `validateNoMixedCaseInput`/`validatePrecalledVcfCoverage`/`extractPrecalledVcfs`/`hasPrecalledSnvVcf`/`hasPrecalledSvVcf`/`hasPrecalledMtVcf`, and a new `test_vcf` profile/pipeline-level test confirming precalled SNV/SV/MT VCFs correctly skip calling while annotation/ranking still run, plus a negative test for the mixed precalled/raw-input case [issue #261](https://github.com/nf-core/raredisease/issues/261)
  [PR #935](https://github.com/nf-core/raredisease/pull/935)
- Add CRAM file input support: accept `cram`/`crai` columns in the samplesheet, converting CRAM to BAM early in the align subworkflow so all downstream tools remain unchanged [issue #261](https://github.com/nf-core/raredisease/issues/261) [PR #933](https://github.com/nf-core/raredisease/pull/933)
- Add `vcf`/`tbi`/`type` columns to the samplesheet schema, enforcing exactly one data type per row (`fastq`, `spring`, `bam`, `cram`, or `vcf`) via a schema-level `oneOf` constraint, and extend `PIPELINE_INITIALISATION` to recognise precalled SNV/SV/MT VCFs and collect them into a new `precalled_vcfs` channel (Only schema/plumbing changes in this PR) [Issue #261](https://github.com/nf-core/raredisease/issues/261) and [PR #927](https://github.com/nf-core/raredisease/pull/927)
- Add `changelog.yml` GitHub Actions workflow to enforce CHANGELOG updates on every PR; PRs can be exempted with the `skip-changelog` label [issue #796](https://github.com/nf-core/raredisease/issues/796) [PR #920](https://github.com/nf-core/raredisease/pull/920)
- Update saltshaker classification reporting by adding customer ID to samples' reports and displaying them as tabs in html [#856](https://github.com/nf-core/raredisease/pull/856)
- Added non-stub tests for `annotate_mt_snvs` [#890](https://github.com/nf-core/raredisease/pull/890)
- Added GATK contamination check for WES/WGS samples as complement to VerifyBamID2, enabled by providing `contamination_sites` and skippable via `--skip_tools gatkcontamination` [#758](https://github.com/nf-core/raredisease/pull/758)
- GATK Contamination results displayed in MultiQC with color-coded thresholds [#758](https://github.com/nf-core/raredisease/pull/758)
- Added non-stub tests for `annotate_mobile_elements` [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #923](https://github.com/nf-core/raredisease/pull/923)
- Added non-stub tests for `call_mobile_elements` [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #924](https://github.com/nf-core/raredisease/pull/924)
- Added non-stub tests for `rank_variants` and `filter_annotate_rank` [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #988](https://github.com/nf-core/raredisease/pull/988)
- Added non-stub tests for `postprocess_MT_calls` and `call_mt_snvs` [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #989](https://github.com/nf-core/raredisease/pull/989)

### `Removed`

- Remove the `test_bam` profile: its all-BAM trio was a strict subset of `test_align`'s coverage, which already exercises both plain-BAM ingestion and CRAM conversion in one run [issue #869](https://github.com/nf-core/raredisease/issues/869) [PR #977](https://github.com/nf-core/raredisease/pull/977)
- Removed the `rtgtools`/`vcfeval` variant-evaluation feature entirely: the `VARIANT_EVALUATION` subworkflow, `rtgtools/format` and `rtgtools/vcfeval` modules, and the `--run_rtgvcfeval`, `--rtg_truthvcfs`, and `--sdf` parameters [issue #963](https://github.com/nf-core/raredisease/issues/963) [PR #964](https://github.com/nf-core/raredisease/pull/964)

### `Changed`

- Add a `tests/lib/TestData.groovy` helper (`TestData.sample('ACC13778A2')`) and use it across the subworkflow nf-tests, replacing 124 repeated inline sample-meta literals [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #1006](https://github.com/nf-core/raredisease/pull/1006)
- Replace the repeated `setup { run("GET_CHROM_SIZES") {…} }` block in eight subworkflow nf-tests (`annotate_genome_snvs`, `annotate_rhocallviz`, `call_mt_snvs`, `call_snv`, `call_snv_deepvariant`, `call_snv_sentieon`, `call_sv_MT`, `postprocess_MT_calls`) with the pre-generated `subworkflow_fixtures/minimal_reference_chrom.sizes` fixture [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #1007](https://github.com/nf-core/raredisease/pull/1007)
- Add `alignment`/`reference_sliced`/`subworkflow_fixtures`/`resources_remapped`/`precalled` path-prefix params to `tests/nextflow.config` and use them across the subworkflow nf-tests, replacing ~560 repeated `params.pipelines_testdata_base_path + '<subdir>/…'` references [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #1005](https://github.com/nf-core/raredisease/pull/1005)
- Migrate `call_snv_MT`'s real and stub tests to the minimal dataset [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #1004](https://github.com/nf-core/raredisease/pull/1004)
- Migrate the stub tests for `call_sv`, `call_sv_manta`, `generate_cytosure_files`, `call_sv_germlinecnvcaller`, and `gens` to the minimal dataset [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #1002](https://github.com/nf-core/raredisease/pull/1002)
- Migrate `contamination`'s real tests to the minimal dataset [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #998](https://github.com/nf-core/raredisease/pull/998)
- Migrate `qc_bam`'s real tests to the minimal dataset [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #999](https://github.com/nf-core/raredisease/pull/999)
- Migrate `prepare_references`'s and `scatter_genome`'s real tests to the minimal dataset [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #996](https://github.com/nf-core/raredisease/pull/996)
- Migrate `call_repeat_expansions`'s and `call_snv`'s real tests to the minimal dataset [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #997](https://github.com/nf-core/raredisease/pull/997)
- Replace the per-profile BWA index memory override (duplicated across every `test_*.config` with a colon-anchored selector) with a single safe floor in `conf/modules/prepare_references.config`, applied unconditionally without affecting `PREPARE_REFERENCES`'s own real-fasta test [issue #869](https://github.com/nf-core/raredisease/issues/869) [PR #975](https://github.com/nf-core/raredisease/pull/975)
- Migrate `call_sv_MT`'s real tests to the minimal dataset [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #995](https://github.com/nf-core/raredisease/pull/995)
- Migrate `call_sv_cnvnator`'s and `call_sv_tiddit`'s real tests to the minimal dataset [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #994](https://github.com/nf-core/raredisease/pull/994)
- Migrate `subsample_mt_frac`'s and `subsample_mt_reads`'s real tests to the minimal dataset [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #993](https://github.com/nf-core/raredisease/pull/993)
- Migrate `call_mobile_elements`'s and `annotate_mobile_elements`'s real tests to the minimal dataset, switching the latter from GRCh37 VEP cache-mode to GTF-mode annotation [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #992](https://github.com/nf-core/raredisease/pull/992)
- Migrate `align_MT`'s real tests to the minimal dataset [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #991](https://github.com/nf-core/raredisease/pull/991)
- Migrate `align_mitochondria`'s and `convert_mt_bam_to_fastq`'s real tests to the minimal dataset [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #990](https://github.com/nf-core/raredisease/pull/990)
- Migrate the `test_singleton` profile's fixture to the minimal 9-region GIAB dataset (same trio, proband alone), matching the default profile's migration; also fixes a stale `skip_tools` override that was silently dropping `smncopynumbercaller` from the test [issue #869](https://github.com/nf-core/raredisease/issues/869) [PR #976](https://github.com/nf-core/raredisease/pull/976)
- Migrate the `test_align` profile's fixture to the minimal dataset (BAM/CRAM trio) [issue #869](https://github.com/nf-core/raredisease/issues/869) [PR #977](https://github.com/nf-core/raredisease/pull/977)
- Migrate the `test_vcf` profile's fixture to the minimal dataset; fixes peddy crashing against the custom sliced-coordinate reference (missing `peddy_sites`) and `genmod models` crashing on the precalled-VCF proband-only samplesheet convention when parent IDs are referenced without a row of their own [issue #869](https://github.com/nf-core/raredisease/issues/869) [PR #978](https://github.com/nf-core/raredisease/pull/978)
- Migrate the `test_sentieon` profile's fixture to the minimal dataset (same trio fastq input as the default profile) [issue #869](https://github.com/nf-core/raredisease/issues/869) [PR #979](https://github.com/nf-core/raredisease/pull/979)
- Replace the default `test` profile's fixture: the chr20/chrM GRCh37 cache-mode dataset is replaced with a 9-region GIAB trio (SNV/SV/STR/ME/MT) sliced to local coordinates on GRCh38, with VEP running in `--gtf` mode. `fastp` and `peddy` are no longer skipped by default; `smncopynumbercaller` is newly skipped (hardcodes SMN1/SMN2 on chr5, not one of this dataset's loci) [issue #869](https://github.com/nf-core/raredisease/issues/869) [PR #972](https://github.com/nf-core/raredisease/pull/972)
- Feed `UPD_SITES`, `UPD_REGIONS`, and `ANNOTATE_RHOCALLVIZ` the unfiltered vcfanno-annotated VCF instead of the VEP/`pre_vep_snv_filter_expression`-filtered one, so stricter filtering (#929) no longer removes variants these subworkflows need [issue #930](https://github.com/nf-core/raredisease/issues/930) [PR #962](https://github.com/nf-core/raredisease/pull/962)
- Bump VEP to 116.1 and default `vep_cache_version` to 116; the updated `ensemblvep/vep` module now takes the VEP cache as a `[meta, path]` tuple instead of a bare path, so `PREPARE_REFERENCES` was updated to emit it in that shape [issue #872](https://github.com/nf-core/raredisease/issues/872) [PR #968](https://github.com/nf-core/raredisease/pull/968)
- Extend the VCF entry point to a fifth type, `repeat` [issue #261](https://github.com/nf-core/raredisease/issues/261) [PR #957](https://github.com/nf-core/raredisease/pull/957)
- Extend the VCF entry point to a fourth type, `me` [issue #261](https://github.com/nf-core/raredisease/issues/261) [PR #955](https://github.com/nf-core/raredisease/pull/955)
- Split `skip_mt_calling` into independently-gated `skip_mt_snv_calling` (gates `CALL_MT_SNVS` only, same behavior as before) and `skip_mt_sv_calling` (new tag gating `CALL_SV_MT` - MitoSalt/SaltShaker and the mitodel/MT-deletion script) [issue #950](https://github.com/nf-core/raredisease/issues/950) [PR #954](https://github.com/nf-core/raredisease/pull/954)
- Extract the clinical-set-filter → CSQ/PLI-annotate → (proband-filter) → rank sequence, previously duplicated once each for SNV/MT/SV/ME in `raredisease.nf`, into a new `FILTER_ANNOTATE_RANK` subworkflow called four times under aliases [issue #952](https://github.com/nf-core/raredisease/issues/952) [PR #953](https://github.com/nf-core/raredisease/pull/953)
- Split mitochondrial SV calling out of `CALL_STRUCTURAL_VARIANTS` into its own top-level `CALL_SV` (nuclear-only) and `CALL_SV_MT` (unchanged) subworkflows called independently from `raredisease.nf`, with the final nuclear+MT SVDB merge moved to the top level; matches how `CALL_SNV`/`CALL_MT_SNVS` are already split; behavior-preserving refactor [issue #948](https://github.com/nf-core/raredisease/issues/948) [PR #951](https://github.com/nf-core/raredisease/pull/951)
- Replace the separate `--sample_id_map` file with an optional `customer_id` column on the main input samplesheet: sample-level customer/external IDs used to label Saltshaker HTML reports and rename VCF2CYTOSURE output are now supplied directly alongside each sample instead of in a second CSV [issue #861](https://github.com/nf-core/raredisease/issues/861) [PR #947](https://github.com/nf-core/raredisease/pull/947)
- Split mitochondrial alignment out of `ALIGN` into its own independently-gated `ALIGN_MITOCHONDRIA` subworkflow, called directly from `raredisease.nf`, matching how `CALL_MT_SNVS` was already split out of `CALL_SNV` [issue #944](https://github.com/nf-core/raredisease/issues/944) [PR #945](https://github.com/nf-core/raredisease/pull/945)
- Migrate all `.set { ch }` and `.tap { ch }` operators to direct channel assignment (`ch = ...`) across `main.nf`, `workflows/raredisease.nf`, and all local subworkflows; behavior-preserving refactor, prerequisite for adopting Nextflow's static type checking [issue #940](https://github.com/nf-core/raredisease/issues/940) [PR #941](https://github.com/nf-core/raredisease/pull/941)
- Speed up `align` subworkflow tests: switch `subworkflows/local/align/tests/main.nf.test` from full GRCh37-scale reference/read data to the small sarscov2 fixtures [issue #938](https://github.com/nf-core/raredisease/issues/938) [PR #939](https://github.com/nf-core/raredisease/pull/939)
- Document the VCF entry point: expand `docs/usage.md`'s samplesheet section with a worked trio example and explicit "what's possible"/"what's not possible" guidance (mixed-input rule, per-type coverage requirement, conflicting-VCF rejection), and add `> **NB**` notes to the SNV/SV/mitochondrial calling sections of `docs/output.md` explaining when those output directories are skipped in favor of the supplied VCF [issue #261](https://github.com/nf-core/raredisease/issues/261)
- Wire up the VCF entry point: a precalled SNV/SV/MT VCF supplied in the samplesheet now auto-skips calling for that type and is substituted directly into annotation. Adds `mt_calling` as a new independently-skippable subworkflow, un-nests `ANNOTATE_STRUCTURAL_VARIANTS` (and its downstream chain) from `skip_sv_calling` to match the existing SNV precedent, fixes `GENERATE_CYTOSURE_FILES` and the `CONCAT_NUCLEAR_AND_MT_SNVS` gate to work correctly with precalled data, and rejects samplesheets where a case mixes precalled VCF rows with fastq/spring/bam/cram rows (a case must be fully precalled or fully processed from raw/aligned reads, not both) [issue #261](https://github.com/nf-core/raredisease/issues/261)[PR #931](https://github.com/nf-core/raredisease/pull/931)
- Split mitochondrial SNV calling out of `CALL_SNV` into its own independently-gated `CALL_MT_SNVS` subworkflow, matching how SV/mobile-element/repeat-expansion calling are each already independently gated in `raredisease.nf`; behavior-preserving refactor, groundwork for letting nuclear and mitochondrial SNV calling be skipped independently [issue #261](https://github.com/nf-core/raredisease/issues/261) [PR #926](https://github.com/nf-core/raredisease/pull/926)
- Un-nest `ANNOTATE_GENOME_SNVS` (and its downstream `GENERATE_CLINICAL_SET_SNV`/`ANN_CSQ_PLI_SNV`/`RANK_VARIANTS_SNV` chain) from `if (!skip_snv_calling)` so nuclear SNV annotation can run independently of whether nuclear SNV calling ran; fix `PEDDY` and `VARIANT_EVALUATION` to route through the unified `ch_call_snv_genome_*` channels instead of referencing `CALL_SNV.out.*` directly (previously unsafe if `snv_calling` were ever skipped), and gate `GENS` on `!skip_snv_calling` since no precalled substitute exists for a gVCF; behavior-preserving refactor, groundwork for the VCF entry point [issue #261](https://github.com/nf-core/raredisease/issues/261) [PR #928](https://github.com/nf-core/raredisease/pull/928)
- Added `--qc_metrics_tool` parameter (`picard` | `riker`) to choose between Picard (`CollectMultipleMetrics` / `CollectHsMetrics` / `CollectWgsMetrics`) and [Riker](https://github.com/fulcrumgenomics/riker) `multi` for BAM QC metrics collection; Riker outputs forwarded to MultiQC [issue #871](https://github.com/nf-core/raredisease/issues/871) [PR #917](https://github.com/nf-core/raredisease/pull/917)
- Consolidate contamination detection into a unified `CONTAMINATION` subworkflow: move `VERIFYBAMID_VERIFYBAMID2` out of `QC_BAM` and `PARSE_CONTAMINATION` out of `raredisease.nf`; replace single `skip_contamination` flag with independent `skip_gatkcontamination` and `skip_verifybamid` flags; `verifybamid` is now skippable via `--skip_tools verifybamid` and auto-skips when `--verifybamid_svd_bed` is not provided; output paths renamed to `contamination/verifybamid/` and `contamination/gatk/` [issue #921](https://github.com/nf-core/raredisease/issues/921) [PR #922](https://github.com/nf-core/raredisease/pull/922)
- Document the `ar x models.bundle` workaround for Sentieon DNAscope `.bundle` model format in `docs/usage.md` and `nextflow_schema.json` [issue #568](https://github.com/nf-core/raredisease/issues/568) [PR #912](https://github.com/nf-core/raredisease/pull/912)
- Pre-resolve the MT analysis flag (`val_analysis_type.matches("wgs|mito") || val_run_mt_for_wes`) into a single named boolean `val_run_mt` in `NFCORE_RAREDISEASE`, replacing `val_run_mt_for_wes` across downstream subworkflow signatures and simplifying repeated conditionals in `PREPARE_REFERENCES`, `ALIGN`, `CALL_SNV`, and `CALL_STRUCTURAL_VARIANTS` [#906](https://github.com/nf-core/raredisease/pull/906)
- Clarify in `docs/usage.md` that the pLI VEP plugin is mandatory when annotation is enabled and LoFtool is optional at the pipeline level [#911](https://github.com/nf-core/raredisease/pull/911)
- Add missing `docs/output.md` sections for GATK contamination check (`qc/contamination/`) and pedigree file (`pedigree/`); fix missing `</details>` closing tag in Peddy section [#904](https://github.com/nf-core/raredisease/pull/904)
- Replace `ch_publish`/`subworkflow_results` with named typed channel emits for fastqc, smncopynumbercaller, peddy, multiqc, and pedigree; remove legacy `publish` emit from `raredisease.nf` [#903](https://github.com/nf-core/raredisease/pull/903)
- Replace `ch_publish`/`subworkflow_results` with named typed channel emits for `annotate_structural_variants` subworkflow [#902](https://github.com/nf-core/raredisease/pull/902)
- Replace `ch_publish`/`subworkflow_results` with named typed channel emits for gens and generate_cytosure_files subworkflows [#899](https://github.com/nf-core/raredisease/pull/899)
- Replace `ch_publish`/`subworkflow_results` with named typed channel emits for `rank_variants` subworkflow [#896](https://github.com/nf-core/raredisease/pull/896)
- Replace `ch_publish`/`subworkflow_results` with named typed channel emits for `variant_evaluation` subworkflow [#897](https://github.com/nf-core/raredisease/pull/897)
- Replace `ch_publish`/`subworkflow_results` with named typed channel emits for `prepare_references` subworkflow [#900](https://github.com/nf-core/raredisease/pull/900)
- Replace `ch_publish`/`subworkflow_results` with named typed channel emits for call_repeat_expansions subworkflow [#893](https://github.com/nf-core/raredisease/pull/893)
- Replace `ch_publish`/`subworkflow_results` with named typed channel emits for call_mobile_elements and annotate_consequence_pli subworkflows; remove `val_publish_dir` parameter from annotate_consequence_pli [#894](https://github.com/nf-core/raredisease/pull/894)
- Replace `ch_publish`/`subworkflow_results` with named typed channel emits for `annotate_mt_snvs` subworkflow [#895](https://github.com/nf-core/raredisease/pull/895)
- Replace `ch_publish`/`subworkflow_results` with named typed channel emits for call_snv, call_snv_deepvariant, and postprocess_MT_calls subworkflows [#863](https://github.com/nf-core/raredisease/pull/863)
- Replace `ch_publish`/`subworkflow_results` with named typed channel emits for annotate_rhocallviz and annotate_genome_snvs subworkflows [#858](https://github.com/nf-core/raredisease/pull/858)
- Expand annotate_rhocallviz test with snapshot assertions [#858](https://github.com/nf-core/raredisease/pull/858)
- Refactor scatter_genome subworkflow: alias GAWK as `GENOME_FAI_TO_BED`, remove `val_save_reference` parameter, move interval flattening into `annotate_genome_snvs` [#857](https://github.com/nf-core/raredisease/pull/857)
- Replace `ch_publish`/`subworkflow_results` with named typed channel emits for qc_bam subworkflow [#853](https://github.com/nf-core/raredisease/pull/853)
- Replace `ch_publish`/`subworkflow_results` with named typed channel emits for alignment and subsample-MT subworkflows [#850](https://github.com/nf-core/raredisease/pull/850)
- Update saltshaker modules to version 1.1.1 so they can run on empty mitosalt output [#856](https://github.com/nf-core/raredisease/pull/856)
- Update metromap to reflect the addition of mitosalt + saltshaker and removal of eklipse [#892](https://github.com/nf-core/raredisease/pull/892)
- Changed default glnexus config from `DeepVariant_unfiltered` to a custom config in `assets/` due to unfixed [bug](https://github.com/dnanexus-rnd/GLnexus/issues/286) [issue #960](https://github.com/nf-core/raredisease/issues/960) [PR #961](https://github.com/nf-core/raredisease/pull/961)
- Change input file for mitosalt from genome-wide fastq to mitochondrial fastq [PR #967](https://github.com/nf-core/raredisease/pull/967)
- Updated tiddit/cov and tiddit/sv to v3.9.7 [PR #1001](https://github.com/nf-core/raredisease/pull/1001)
- Updated vcf2cytosure to v0.10.0 [PR #1003](https://github.com/nf-core/raredisease/pull/1003)
- Updated `deepvariant/rundeepvariant` to v1.10.0 [PR #1010](https://github.com/nf-core/raredisease/pull/1010)

### `Fixed`

- Fix `call_sv`'s standalone subworkflow test failing against the new minimal dataset with "Minimum memory limit allowed is 6MB": `BWA_INDEX`'s default memory (proportional to fasta size) computes below Docker's floor for the tiny sliced reference; the test's own setup step now applies the same `[6.B * fasta.size(), 100.MB].max()` floor already used for the real `PREPARE_REFERENCES:BWA_INDEX_GENOME` invocation in `conf/modules/prepare_references.config` [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #983](https://github.com/nf-core/raredisease/pull/983)
- Fix `call_snv_sentieon`'s standalone test passing the genome fasta and fai in the wrong argument order (a pre-existing bug masked by stub mode never touching file content) [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #987](https://github.com/nf-core/raredisease/pull/987)
- Fix `call_snv_deepvariant`'s standalone test hardcoding `--regions="chr22:0-40001"`, left over from the old dataset; the new minimal dataset has no chr22 region, so DeepVariant errored with "regions to call is empty" [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #987](https://github.com/nf-core/raredisease/pull/987)
- Fix `call_snv_deepvariant`'s standalone test config missing `-c CHROM,FROM,TO,FOUND_IN` on `BCFTOOLS_ANNOTATE` (present in the real pipeline config but never mirrored into the subworkflow test's own `tests/nextflow.config`), which made the real test's `bcftools annotate` step fail with "The -c option not given" [issue #795](https://github.com/nf-core/raredisease/issues/795) [PR #987](https://github.com/nf-core/raredisease/pull/987)
- Fix pipeline-level `conf/modules/*.config` selectors leaking into subworkflow-level nf-test isolation: `withName` patterns nested under a locally-tested subworkflow (e.g. `.*CALL_SV:CALL_SV_MANTA:MANTA`) matched that subworkflow's own standalone test run as well as the real pipeline, since `.*` can match an empty prefix; anchored ~165 selectors across 30 `conf/modules/*.config` files to `.*:PARENT:CHILD:PROC` so they only apply when reached via a genuine outer (pipeline) prefix, and backfilled the subworkflow-local `tests/nextflow.config` overrides (adding the required `config "./nextflow.config"` test directive where missing) that several subworkflows had been unknowingly relying on the leaked pipeline values for [issue #981](https://github.com/nf-core/raredisease/issues/981) [PR #982](https://github.com/nf-core/raredisease/pull/982)
- Fix `pre_vep_snv_filter_expression`'s default hardcoding `GNOMADAF_popmax`: gnomAD v4 renamed the field to `GNOMADAF_grpmax`, so this filter failed with "tag not defined in VCF header" against any gnomAD v4-labeled annotation, not a dataset-specific issue [issue #869](https://github.com/nf-core/raredisease/issues/869) [PR #969](https://github.com/nf-core/raredisease/pull/969)
- Fix `ANNOTATE_MOBILE_ELEMENTS:BCFTOOLS_VIEW_FILTER` discarding every mobile-element call: RetroSeq never sets `FILTER=PASS` (always `.`), so `--apply-filters PASS` kept nothing; now accepts `.` as well as `PASS` [issue #869](https://github.com/nf-core/raredisease/issues/869) [PR #969](https://github.com/nf-core/raredisease/pull/969)
- Fix `ENSEMBLVEP_SNV`/`ENSEMBLVEP_MT`'s SpliceAI plugin argument scoring every indel against SNV-only data: `indel=` was copy-pasted from `snv=` and pointed at the SNV-scores file for both. Also fixes the same copy-paste bug in `annotate_mt_snvs`'s standalone subworkflow-test config, whose real (non-stub) test never staged an indel-named SpliceAI file either, and unifies the SpliceAI plugin filenames used across cache mode and `--gtf` mode onto `spliceai_snv.vcf.gz`/`spliceai_indel.vcf.gz`, backed by [nf-core/test-datasets#2212](https://github.com/nf-core/test-datasets/pull/2212) renaming the shared `reference/` files [issue #869](https://github.com/nf-core/raredisease/issues/869) [PR #969](https://github.com/nf-core/raredisease/pull/969)
- Fix `MERGE_NUCLEAR_AND_MT_SVS` crashing with "If priority is used, one tag per VCF is needed" whenever both nuclear and mitochondrial SV calling ran: `CALL_SV` no longer pre-merges the nuclear caller VCFs (tiddit/manta/gcnvcaller/cnvnator) into one file before the top-level merge; it now emits the individual per-caller VCFs so the single merge in `raredisease.nf` always has one priority tag per input VCF. Regression from the `CALL_SV`/`CALL_SV_MT` split in [issue #948](https://github.com/nf-core/raredisease/issues/948) [PR #959](https://github.com/nf-core/raredisease/pull/959)
- Speed up and stabilise the call_snv_deepvariant/call_snv (DeepVariant) subworkflow tests by using small chr22 fixtures, adding 2-sample coverage for GLnexus joint genotyping, and restoring the deterministic variantsMD5 snapshot assertion. [#918](https://github.com/nf-core/raredisease/pull/918)
- Emit an error at startup when `vep_filters_scout_fmt` or `vep_filters` contains no records (headers or empty lines only), which would otherwise cause the clinical set to silently contain 0 variants [#913](https://github.com/nf-core/raredisease/pull/913)
- Add missing CADD 1.7.3 module update to the v3.0.0 `Tool updates` table in `CHANGELOG.md` [issue #888](https://github.com/nf-core/raredisease/issues/888) [PR #919](https://github.com/nf-core/raredisease/pull/919)
- Add changelog entry requirement (including the Parameters table) to the contribution workflow in `CONTRIBUTING.md` [issue #797](https://github.com/nf-core/raredisease/issues/797) [PR #919](https://github.com/nf-core/raredisease/pull/919)
- Fix `--hisat2` parameter being declared but never consumed, causing the HISAT2 genome index to always be rebuilt [#905](https://github.com/nf-core/raredisease/pull/905)
- Fix inconsistent sample column order in Sentieon SNV family VCF by sorting per-sample VCFs by filename before merging, consistent with DeepVariant and MT paths [#908](https://github.com/nf-core/raredisease/pull/908)
- Fix intermittent `CALL_SNV_DEEPVARIANT - wgs` test failure caused by non-deterministic GLnexus quality scores by replacing `variantsMD5` with `vcf.summary` [#850](https://github.com/nf-core/raredisease/pull/850)
- Fix swapped `run_mt_for_wes`/`skip_split_multiallelics` arguments in the `CALL_SNV` call, which disabled multiallelic splitting (and inverted MT-for-WES) when `--run_mt_for_wes` was set [#854](https://github.com/nf-core/raredisease/issues/854)

### Parameters

| Old parameter | New parameter                 |
| ------------- | ----------------------------- |
|               | contamination_sites           |
|               | contamination_sites_tbi       |
|               | pre_vep_snv_filter_expression |
|               | glnexus_config                |
|               | vep_gtf                       |
|               | vep_gtf_tbi                   |
|               | peddy_sites                   |
|               | duplicates_marker             |
|               | somalier_sites_vcf            |

### Tool updates

| Tool                         | Old version | New version |
| ---------------------------- | ----------- | ----------- |
| gatk4/calculatecontamination |             | 4.6.2.0     |
| gatk4/getpileupsummaries     |             | 4.6.2.0     |
| Saltshaker                   | 1.0.0       | 1.1.1       |
| Ensemblvep                   | 110.1       | 116.1       |
| tiddit/sv                    | 3.9.5       | 3.9.7       |
| tiddit/cov                   | 3.9.5       | 3.9.7       |
| vcf2cytosure                 | 0.9.3       | 0.10.0      |
| deepvariant                  | 1.9.0       | 1.10.0      |

## 3.1.2 - Princess Peach (patch) [2026-07-06]

### `Fixed`

- Fix `svdb/merge` mislabelling VCF caller tags by passing `sort_inputs = false` to `SVDB_MERGE` in `call_structural_variants`; the VCF list is already in the correct caller-priority order from the upstream `concat` chain, so in-module re-sorting by filename was causing tag misassignment [#910](https://github.com/nf-core/raredisease/pull/910)
- Changed filter logic for merging mitochondrial snvs to report variants which pass filter in at least one sample [#915](https://github.com/nf-core/raredisease/pull/915)

## 3.1.1 - Princess Peach (patch) [2026-06-24]

### `Fixed`

- Patch `deepvariant/rundeepvariant` to tee stdout/stderr to a log file and exit non-zero when `queue.Empty` or `BrokenPipeError` is detected, catching silent failures that previously caused the process to appear successful [#889](https://github.com/nf-core/raredisease/pull/889)

## 3.1.0 - Princess Peach [2026-06-16]

### `Added`

- Parameter `cadd_prescored` to pass a directory of pre-scored CADD indel annotations to the CADD process in genome and mitochondrial SNV annotation subworkflows [#866](https://github.com/nf-core/raredisease/pull/866)
- Parameter `manta_call_regions` to restrict Manta SV calling to specified regions (e.g. primary chromosomes) via a bgzipped, tabix-indexed BED file, reducing runtime without affecting other callers [#867](https://github.com/nf-core/raredisease/pull/867)
- Local `FILTERVEP` module using a Python reimplementation of Ensembl's `filter_vep`, replacing the `ENSEMBLVEP_FILTERVEP` module with a lighter cyvcf2-based alternative [#870](https://github.com/nf-core/raredisease/pull/870)
- `bwafastalign/index` nf-core module and `bwafastalign` parameter to support index preparation for the bwa-fastalign genome aligner [#877](https://github.com/nf-core/raredisease/pull/877)
- `bwafastalign/mem` nf-core module to support genome alignment with bwa-fastalign when `--aligner bwafastalign` is set [#880](https://github.com/nf-core/raredisease/pull/880)

### `Changed`

- Replace `ENSEMBLVEP_FILTERVEP` with local `FILTERVEP` in the clinical set subworkflow, renamed from `VCF_FILTER_BCFTOOLS_ENSEMBLVEP` to `VCF_FILTER_BCFTOOLS_FILTERVEP` [#870](https://github.com/nf-core/raredisease/pull/870)
- Increase default mbuffer memory value from 3GB to 8GB [#880](https://github.com/nf-core/raredisease/pull/880)
- Update `bwameme/mem` to new nf-core module signature: `val mbuffer` and `val samtools_threads` replaced by `ext.args2` and `ext.args3` [#881](https://github.com/nf-core/raredisease/pull/881)

### `Fixed`

- Add a bcftools norm split-multiallelics step after merging standard and shifted MT calls to handle new multiallelic sites introduced by bcftools merge [#855](https://github.com/nf-core/raredisease/pull/855)

### Parameters

| Old parameter | New parameter          |
| ------------- | ---------------------- |
|               | bwafastalign           |
|               | cadd_prescored         |
|               | manta_call_regions     |
|               | manta_call_regions_tbi |

### Tool updates

| Tool          | Old version | New version |
| ------------- | ----------- | ----------- |
| bwa-fastalign |             | 1.0.0       |
| saltshaker    | 1.0.0       | 1.1.1       |

## 3.0.0 - Mario [2026-05-12]

### `Added`

- Interval parameter in the default retroseq call [#717](https://github.com/nf-core/raredisease/pull/717)
- Tests for call_repeat_expansions and qc_bam subworkflows [#713](https://github.com/nf-core/raredisease/pull/713)
- Feature to subsample mitochondrial alignments based on number of reads [#748](https://github.com/nf-core/raredisease/pull/748)
- Functionality to generate coverage information using Sambamba depth [#752](https://github.com/nf-core/raredisease/pull/752)
- Parameter to pass a file containing new sample ids to use with multiqc [#764](https://github.com/nf-core/raredisease/pull/764)
- A helper function channelFromPath to create channels in a readable fashion in main.nf [#766](https://github.com/nf-core/raredisease/pull/766)
- A helper function channelFromPathWithMeta to create channels in a readable fashion in main.nf [#767](https://github.com/nf-core/raredisease/pull/767)
- A helper function channelFromSamplesheet to create channels in a readable fashion in main.nf [#767](https://github.com/nf-core/raredisease/pull/767)
- A parameter `homoplasmy_af_threshold` to set genotypes of MT SNVs to 1/1 (homoplasmic) when AF >=`homoplasmy_af_threshold` [#768](https://github.com/nf-core/raredisease/pull/768)
- Topic channels to local modules to caputure versions [#774](https://github.com/nf-core/raredisease/pull/774)
- MitoSalt to detect mitochondrial deletions [#743](https://github.com/nf-core/raredisease/pull/743)
- Tests for some of the subworkflows [#780](https://github.com/nf-core/raredisease/pull/780)
- Tests for some of the subworkflows [#782](https://github.com/nf-core/raredisease/pull/782)
- Tests for some of the subworkflows [#783](https://github.com/nf-core/raredisease/pull/783)
- Test tags for dependent modules in subworkflow tests [#800](https://github.com/nf-core/raredisease/pull/800)
- Stub test for scatter_genome subworkflow [#802](https://github.com/nf-core/raredisease/pull/802)
- Add CAT_FASTQ before SEQTK_SAMPLE in call_sv_MT to merge reads across lanes before subsampling for MitoSalt [#799](https://github.com/nf-core/raredisease/pull/799)
- Add new local SPLIT_CHR module to split reference FASTA by chromosome for CNVnator, and pass genome flag to all CNVnator steps [#799](https://github.com/nf-core/raredisease/pull/799)
- Add peddy --sites hg38 argument when running with GRCh38 [#799](https://github.com/nf-core/raredisease/pull/799)
- Saltshaker for downstream processing of mitochondrial SV calls from MitoSAlt [#775](https://github.com/nf-core/raredisease/pull/775)
- Env variable NXF_SINGULARITY_NEW_PID_NAMESPACE = false to accommodate hisat2 running with latest Nextflow and Singularity [#775](https://github.com/nf-core/raredisease/pull/775)
- Parameter `exclude_alt` to filter alignments to alt/unplaced contigs after alignment using samtools view, retaining only primary chromosomes (GRCh37: 1-22,X,Y,MT / GRCh38: chr1-chr22,chrX,chrY,chrM). Note that enabling this will restrict variant calling to these chromosomes [#803](https://github.com/nf-core/raredisease/pull/803)
- Stub test for all the remaning subworkflows that were lacking it: align_bwa_bwamem2_bwameme, align_MT, align (bwameme - wes), align_sentieon, call_repeat_expansions, prepare_references, qc_bam [#820](https://github.com/nf-core/raredisease/pull/820)
- Parameters `save_all_mapped_as_cram` and `save_noalt_mapped_as_cram` to replace `save_mapped_as_cram`, allowing independent control over publishing unfiltered and alt-filtered alignment files as CRAM [#807](https://github.com/nf-core/raredisease/pull/807)
- Parameter `run_vcfanno_db_sanity_check` to check vcfanno database files for zero records and remove the corresponding annotation blocks from the TOML config before running vcfanno [#821](https://github.com/nf-core/raredisease/pull/821)
- Added `--skip_split_multiallelics` parameter to allow users to skip the `bcftools norm --multiallelics -both` step in SNV calling (DeepVariant and Sentieon), which can cause indel quality degradation in single-interval runs [#823](https://github.com/nf-core/raredisease/pull/823)
- Add find/concatenate step to concatenate saltshaker classification files before creating the html report, so the final report is case-level. [#826](https://github.com/nf-core/raredisease/pull/826)
- Extended vcfanno database sanity check to include extra vcfanno resources (`vcfanno_extra`) alongside the main resources, and moved the check upstream to `raredisease.nf` so it covers both genome and mitochondrial SNV annotation subworkflows [#834](https://github.com/nf-core/raredisease/pull/834)
- Add full test to call_sv_MT subworkflow [#874](https://github.com/nf-core/raredisease/pull/874)

### `Changed`

- Use distinct output filenames for bcfools (in call_mobile_elements subworkflow) and svdb (in call_sv_tiddit subworkflow) [#716](https://github.com/nf-core/raredisease/pull/716)
- Use nf-core's most severe consequence & pli scripts instead of local ones [#732](https://github.com/nf-core/raredisease/pull/732)
- Use nf-core's VCF_FILTER_BCFTOOLS_ENSEMBLVEP subworkflow to generate clinical set instead of a local subworkflow [#727](https://github.com/nf-core/raredisease/pull/727)
- Don't call mobile elements in mitochondrial DNA. [#741](https://github.com/nf-core/raredisease/pull/741)
- Call SVs in mitochondria using mitochondrial alignments in the genome alignment files instead of from BAM files generated by the mitochondrial subworkflow. [#742](https://github.com/nf-core/raredisease/pull/742)
- Update gens-preproc script [#747](https://github.com/nf-core/raredisease/pull/747)
- Removed parameter `bwa_as_fallback` [#763](https://github.com/nf-core/raredisease/pull/763)
- Sambamba depth now filters on not duplicates and not failed_quality_control [#768](https://github.com/nf-core/raredisease/pull/768)
- Removed eKLIPse [#743](https://github.com/nf-core/raredisease/pull/743)
- Removed haplocheck [#778](https://github.com/nf-core/raredisease/pull/778)
- Removed HmtNote [#779](https://github.com/nf-core/raredisease/pull/779)
- Updated svbd module [#781](https://github.com/nf-core/raredisease/pull/781)
- Migrate file publishing from publishDir to a centralized output {} block for some workflows [#784](https://github.com/nf-core/raredisease/pull/784)
- Replace local gens module with nf-core module [#785](https://github.com/nf-core/raredisease/pull/785)
- Migrate file publishing from publishDir to a centralized output {} block for some workflows [#787](https://github.com/nf-core/raredisease/pull/787)
- Migrate file publishing from publishDir to a centralized output {} block for some workflows [#788](https://github.com/nf-core/raredisease/pull/788)
- Migrate file publishing from publishDir to a centralized output {} block for some workflows [#789](https://github.com/nf-core/raredisease/pull/789)
- Remove redundant TABIX processes, and update configs for nf-test [#790](https://github.com/nf-core/raredisease/pull/790)
- Migrate file publishing from publishDir to a centralized output {} block for some workflows [#791](https://github.com/nf-core/raredisease/pull/791)
- Remove redundant ZIP_TABIX steps after VCFANNO in annotate_genome_snvs and annotate_mt_snvs by using VCFANNO's direct tbi output [#799](https://github.com/nf-core/raredisease/pull/799)
- Collect genome fasta/fai channel in call_sv_tiddit to prevent per-sample re-emission [#799](https://github.com/nf-core/raredisease/pull/799)
- Update cadd_resources channel to use channelFromPathWithMeta and set channelFromSamplesheet calls for svdb/ME resources as non-mandatory [#799](https://github.com/nf-core/raredisease/pull/799)
- Run MitoSAlt.pl from bin rather than within container [#775](https://github.com/nf-core/raredisease/pull/775)
- Include mitochonrdial SV calls in combined SV vcf, change call_sv output directory structure to remove mitochondria/ and genome/ [#775](https://github.com/nf-core/raredisease/pull/775)
- Remove Qualimap and Haplogrep3 as they were made redundant by Picard and VerifyBamID2 [#801](https://github.com/nf-core/raredisease/pull/801)
- Remove env variable NXF_SINGULARITY_NEW_PID_NAMESPACE from the config since this has to be set outside the subworkflow [#804](https://github.com/nf-core/raredisease/pull/804)
- Run UPD_SITES, UPD_REGIONS, and CHROMOGRAPH for UPD only when analysis type is WGS [#806](https://github.com/nf-core/raredisease/pull/806)
- Change saltshaker classification output from txt to html [#808](https://github.com/nf-core/raredisease/pull/808)
- Sort parameters of `CALL_STRUCTURAL_VARIANTS` and `CALL_SV_MANTA` alphabetically [#821](https://github.com/nf-core/raredisease/pull/821)

### `Fixed`

- Fixed argument order of `ch_genome_fai` and `ch_genome_fasta` in the `CALL_SNV_SENTIEON` subworkflow [#811](https://github.com/nf-core/raredisease/pull/811)
- Ensure deterministic sample ordering in Manta SV output by sorting BAM/BAI channel inputs [#815](https://github.com/nf-core/raredisease/pull/815)
- Fixed inconsistencies in JSON schema [#714](https://github.com/nf-core/raredisease/pull/714)
- Fixed conda declaration in the add_varcallername_to_bed module [#733](https://github.com/nf-core/raredisease/pull/733)
- Fixed CADD annotation to support chr prefix [#745](https://github.com/nf-core/raredisease/pull/745)
- Fixed mismatch between VCF and ROH calls when analysing multiple samples [#755](https://github.com/nf-core/raredisease/pull/755)
- Fixed pipeline to run variant calling even with bait_padding set to 0 [#757](https://github.com/nf-core/raredisease/pull/757)
- Fixed mitosalt channel handling so it runs on all samples in a trio [#826](https://github.com/nf-core/raredisease/pull/826)
- Fixed runtime errors in `call_sv_MT` and `call_structural_variants` when MitoSAlt produces no structural variant calls [#837](https://github.com/nf-core/raredisease/pull/837)
- Fixed vcfanno sanity check map closure to handle `ch_vcfanno_resources` emitting a list of paths [#837](https://github.com/nf-core/raredisease/pull/837)
- Fixed `PREP_MITOSALT` msconfig output being consumed as a queue channel by converting it to a value channel with `.collect()` before passing to `MITOSALT` [#837](https://github.com/nf-core/raredisease/pull/843)

### Parameters

| Old parameter       | New parameter               |
| ------------------- | --------------------------- |
|                     | sambamba_regions            |
| bwa_as_fallback     |                             |
|                     | multiqc_samples             |
|                     | homoplasmy_af_threshold     |
|                     | exclude_alt                 |
| save_mapped_as_cram |                             |
|                     | save_all_mapped_as_cram     |
|                     | save_noalt_mapped_as_cram   |
|                     | run_vcfanno_db_sanity_check |
|                     | skip_split_multiallelics    |

### Tool updates

| Tool                  | Old version | New version |
| --------------------- | ----------- | ----------- |
| bcftools              | 1.20        | 1.21        |
| bwa                   | 0.7.18      | 0.7.19      |
| cadd                  | 1.6.post1   | 1.7.3       |
| deepvariant           | 1.8.0       | 1.9.0       |
| eKLIPse               | 1.8         |             |
| ensemblvep/vep        | 110         | 110.1       |
| ensemblvep/filtervep  | 113         | 115.2       |
| fastp                 | 0.23.4      | 1.0.1       |
| gatk4                 | 4.5.0.0     | 4.6.2.0     |
| gawk                  | 5.3.0       | 5.3.1       |
| genmod                | 3.9         | 3.10.2      |
| gens-preproc          | 1.0.11      |             |
| gens/preparecovandbaf |             | 1.1.5       |
| haplocheck            | 1.3.3       |             |
| haplogrep3            | 3.2.2       |             |
| hmtnote               | 0.7.2       |             |
| htslib                | 1.20        | 1.21        |
| MitoSalt              |             | 1.1.1       |
| mosdepth              | 0.3.8       | 0.3.11      |
| multiqc               | 1.32        | 1.33        |
| ngsbits               | 202411      | 202512      |
| picard                | 3.3.0       | 3.4.0       |
| pigz                  | 2.3.4       | 2.8         |
| qualimap              | 2.3         |             |
| sambamba              |             | 1.0.1       |
| samtools              | 1.21        | 1.22.1      |
| sentieon              | 202503      | 202503.02   |
| stranger              | 0.9.4       | 0.10.0      |
| svdb                  | 2.8.3       | 2.8.4       |
| tiddit                | 3.6.1       | 3.9.5       |
| ucsc                  | 447         | 482         |
| vcfanno               | 0.3.5       | 0.3.7       |
| vcf2cytosure          | 0.9.1       | 0.9.3       |

## 2.6.0 - Cacofonix [2025-06-25]

### `Added`

- A feature to start the workflow from duplicate-marked bam files [#682](https://github.com/nf-core/raredisease/pull/682)
- A functionality to take spring files as input [#678](https://github.com/nf-core/raredisease/pull/678)

### `Changed`

- Refactored code to only handle clinical variants in the generate_clinical_set workflow [#693](https://github.com/nf-core/raredisease/pull/693)
- Refactored `schema_input.json` and `nextflow_schema.json` to improve the error messages and validations of the pipeline [#692](https://github.com/nf-core/raredisease/pull/692)
- Updated `add_most_severe_consequence` and `add_most_severe_pli` to fix spelling and language server warnings [#689](https://github.com/nf-core/raredisease/pull/689)
- Refactored code to address issues highlighted by language server [#688](https://github.com/nf-core/raredisease/pull/688)
- Changed for loop to each in create_pedigree_file [#683](https://github.com/nf-core/raredisease/pull/683)
- Refactored repeat annotation and updated Stranger [#702](https://github.com/nf-core/raredisease/pull/702)

### `Fixed`

- Errors due to channel name and structure inconsistencies in the sentieon SNV calling subworkflow[#688](https://github.com/nf-core/raredisease/pull/688)

### Parameters

| Old parameter            | New parameter     |
| ------------------------ | ----------------- |
| skip_haplogrep3          | skip_tools        |
| skip_fastp               |                   |
| skip_gens                |                   |
| skip_germlinecnvcaller   |                   |
| skip_peddy               |                   |
| skip_smncopynumbercaller |                   |
| skip_vcf2cytosure        |                   |
| skip_me_calling          | skip_subworkflows |
| skip_me_annotation       |                   |
| skip_mt_annotation       |                   |
| skip_mt_subsample        |                   |
| skip_repeat_annotation   |                   |
| skip_repeat_calling      |                   |
| skip_snv_annotation      |                   |
| skip_snv_calling         |                   |
| skip_sv_annotation       |                   |
| skip_sv_calling          |                   |

### Tool updates

| Tool                        | Old version | New version |
| --------------------------- | ----------- | ----------- |
| DeepVariant                 | 1.6.1       | 1.8.0       |
| add_most_severe_consequence | 1.0         | 1.1         |
| add_most_severe_pli         | 1.0         | 1.1         |
| stranger                    | 0.9.2       | 0.9.4       |

## v2.5.0 - Fulliautomatix [2025-05-22]

### `Added`

- A new parameter `concatenate_snv_calls` to generate a concatenated VCF file containing unannotated nuclear & mitochondrial SNV calls [#699](https://github.com/nf-core/raredisease/pull/699)
- Functionality to check contamination in samples using VerifyBamID2 [#701](https://github.com/nf-core/raredisease/pull/701)
- New parameters `verifybamid_svd_bed`, `verifybamid_svd_mu`, and `verifybamid_svd_ud` to supply reference files for VerifyBamID2 [#701](https://github.com/nf-core/raredisease/pull/701)

### `Changed`

- Default to remove mitochondrial variants with FILTER status not equal to PASS [#697](https://github.com/nf-core/raredisease/pull/697)

### `Fixed`

- Sort the input files before vcf2cytosure is invoked [#697](https://github.com/nf-core/raredisease/pull/697)
- Use '--mitochondria-mode' by default when running Gatk4 FilterMutectCalls on mitochondrial variants[#697](https://github.com/nf-core/raredisease/pull/697)

### Parameters

| Old parameter | New parameter         |
| ------------- | --------------------- |
|               | concatenate_snv_calls |
|               | verifybamid_svd_bed   |
|               | verifybamid_svd_mu    |
|               | verifybamid_svd_ud    |

### Tool updates

| Tool         | Old version | New version |
| ------------ | ----------- | ----------- |
| VerifyBamID2 |             | 2.0.1       |

## 2.4.0 - Vitalstatistix [2025-02-24]

### `Added`

- Add markduplicates metrics to multiqc [#679](https://github.com/nf-core/raredisease/pull/679)

### `Changed`

- Update SVDB/merge to fix language server problems [#684](https://github.com/nf-core/raredisease/pull/684)

### `Fixed`

### Parameters

| Old parameter | New parameter |
| ------------- | ------------- |
|               |               |

### Tool updates

| Tool | Old version | New version |
| ---- | ----------- | ----------- |
|      |             |             |

## 2.3.0 - Getafix [2025-02-13]

### `Added`

- A new option `skip_haplogrep3` to skip haplogrep3 [#675](https://github.com/nf-core/raredisease/pull/675)
- A new analysis option `mito` to call and annotate only mitochondrial variants [#608](https://github.com/nf-core/raredisease/pull/608)
- An option `extract_alignments` to restrict analysis to specific contigs [#644](https://github.com/nf-core/raredisease/pull/644)
- Fastp and ngsbits output files as input of MultiQC [#647](https://github.com/nf-core/raredisease/pull/647/).
- Haplocheck output file as input of MultiQC [#662](https://github.com/nf-core/raredisease/pull/662)

### `Changed`

- Only Structural variants that PASS the filter will be reported and annotated [#673](https://github.com/nf-core/raredisease/pull/673)
- Update haplogrep to v3.2.2 [#672](https://github.com/nf-core/raredisease/pull/672)
- d4 files are not generated by default anymore [#648](https://github.com/nf-core/raredisease/pull/648)
- Suffix used to identify unique fastq pairs from "\_T" to "\_LNUMBER" [#638](https://github.com/nf-core/raredisease/pull/638)
- Merge output from germlinecnvcaller [#635](https://github.com/nf-core/raredisease/pull/635)
- Update tools [#623](https://github.com/nf-core/raredisease/pull/623)
- Update output file name prefix for upd and chromograph to sample-based [#620](https://github.com/nf-core/raredisease/pull/620)
- Update tools [#619](https://github.com/nf-core/raredisease/pull/619)
- Report only variants above 5% heteroplasmy in the clinical vcf file for mitochondria [#616](https://github.com/nf-core/raredisease/pull/616)

### `Fixed`

- Download tests [#667](https://github.com/nf-core/raredisease/pull/667)
- Restrict deepvariant analysis of WES samples to bait regions [#633](https://github.com/nf-core/raredisease/pull/633), [#658](https://github.com/nf-core/raredisease/pull/658)
- bcftools annotate declaration in annotate CADD subworkflow [#624](https://github.com/nf-core/raredisease/pull/624)
- Rhocallviz subworkflow will only be invocated once per sample [#621](https://github.com/nf-core/raredisease/pull/621)
- Updated createCaseChannel function to include a check for maternal and paternal ids being set to a numeric 0 [#643](https://github.com/nf-core/raredisease/pull/643)
- Fixed issue of parsing sample sex in the configs [#659](https://github.com/nf-core/raredisease/pull/659)
- Fixed how meta.id was set for input files, giving the resulting files clearer filenames [#661](https://github.com/nf-core/raredisease/pull/661)

### Parameters

| Old parameter | New parameter       |
| ------------- | ------------------- |
|               | extract_alignments  |
|               | restrict_to_contigs |
|               | skip_haplogrep3     |

### Tool updates

| Tool       | Old version | New version |
| ---------- | ----------- | ----------- |
| bcftools   | 1.18        | 1.20        |
| ensemblvep | 112         | 110         |
| genmod     | 3.8.2       | 3.9         |
| haplogrep  | 2.4.0       | 3.2.2       |
| mosdepth   | 0.3.6       | 0.3.8       |
| multiqc    | 1.21        | 1.26        |
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
