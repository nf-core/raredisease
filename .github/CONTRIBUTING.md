# `nf-core/raredisease`: Contributing Guidelines

Hi there!
Many thanks for taking an interest in improving nf-core/raredisease.

We try to manage the required tasks for nf-core/raredisease using GitHub issues, you probably came to this page when creating one.
Please use the pre-filled template to save time.

However, don't be put off by this template - other more general issues and suggestions are welcome!
Contributions to the code are even more welcome ;)

> [!NOTE]
> If you need help using or modifying nf-core/raredisease then the best place to ask is on the nf-core Slack [#raredisease](https://nfcore.slack.com/channels/raredisease) channel ([join our Slack here](https://nf-co.re/join/slack)).

## Table of contents

- [General](#general)
  - [Contribution workflow](#contribution-workflow)
  - [Running tests](#running-tests)
  - [Patch](#patch)
  - [Getting help](#getting-help)
  - [Nextflow version bumping](#nextflow-version-bumping)
  - [Images and figures](#images-and-figures)
  - [GitHub Codespaces](#github-codespaces)
- [Pipeline-specific conventions](#pipeline-specific-conventions)
  - [Architecture & structure](#architecture--structure)
  - [Adding a new step](#adding-a-new-step)
  - [Channel conventions](#channel-conventions)
  - [Params & analysis types](#params--analysis-types)
  - [Publishing](#publishing)
  - [Configuration](#configuration)
  - [Writing tests](#writing-tests)
  - [Style](#style)
  - [Adding citations](#adding-citations)

## General

### Contribution workflow

1. Check that there isn't already an issue about your idea in the [nf-core/raredisease issues](https://github.com/nf-core/raredisease/issues) to avoid duplicating work. If there isn't one already, please create one so that others know you're working on this.
2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [nf-core/raredisease repository](https://github.com/nf-core/raredisease) to your GitHub account.
3. Make the necessary changes / additions within your forked repository following the conventions below.
4. Use `nf-core pipelines schema build` to add any new parameters to `nextflow_schema.json`.
5. Submit a Pull Request against the `dev` branch and wait for the code to be reviewed and merged.

If you're not used to this workflow with git, you can start with some [docs from GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or even their [excellent `git` resources](https://try.github.io/).

### Running tests

You have the option to test your changes locally by running the pipeline. For receiving warnings about process selectors and other `debug` information, it is recommended to use the debug profile. Execute all the tests with the following command:

```bash
nextflow run . -profile debug,test,docker --outdir <OUTDIR>
```

When you create a pull request with changes, [GitHub Actions](https://github.com/features/actions) will run automatic tests.
Typically, pull-requests are only fully reviewed when these tests are passing, though of course we can help out before then.

There are typically two types of tests that run:

#### Lint tests

`nf-core` has a [set of guidelines](https://nf-co.re/developers/guidelines) which all pipelines must adhere to.
To enforce these and ensure that all pipelines stay in sync, we have developed a helper tool which runs checks on the pipeline code. This is in the [nf-core/tools repository](https://github.com/nf-core/tools) and once installed can be run locally with the `nf-core pipelines lint <pipeline-directory>` command.

If any failures or warnings are encountered, please follow the listed URL for more documentation.

#### Pipeline tests

Each `nf-core` pipeline should be set up with a minimal set of test-data.
`GitHub Actions` then runs the pipeline on this data to ensure that it exits successfully.
If there are any failures then the automated tests fail.
These tests are run both with the latest available version of `Nextflow` and also the minimum required version that is stated in the pipeline code.

### Patch

:warning: Only in the unlikely and regretful event of a release happening with a bug.

- On your own fork, make a new branch `patch` based on `upstream/main` or `upstream/master`.
- Fix the bug, and bump version (X.Y.Z+1).
- Open a pull-request from `patch` to `main`/`master` with the changes.

### Getting help

For further information/help, please consult the [nf-core/raredisease documentation](https://nf-co.re/raredisease/usage) and don't hesitate to get in touch on the nf-core Slack [#raredisease](https://nfcore.slack.com/channels/raredisease) channel ([join our Slack here](https://nf-co.re/join/slack)).

### Nextflow version bumping

If you are using a new feature from core Nextflow, bump the minimum required version with:

```bash
nf-core pipelines bump-version --nextflow . [min-nf-version]
```

### Images and figures

For overview images and other documents we follow the nf-core [style guidelines and examples](https://nf-co.re/developers/design_guidelines).

### GitHub Codespaces

This repo includes a devcontainer configuration which will create a GitHub Codespaces for Nextflow development! This is an online developer environment that runs in your browser, complete with VSCode and a terminal.

To get started:

- Open the repo in [Codespaces](https://github.com/nf-core/raredisease/codespaces)
- Tools installed
  - nf-core
  - Nextflow

Devcontainer specs:

- [DevContainer config](.devcontainer/devcontainer.json)

## Pipeline-specific conventions

### Architecture & structure

- **One subworkflow per biological task** — alignment, QC, variant calling, annotation, and ranking are each their own subworkflow under `subworkflows/local/`. Don't add logic to `workflows/raredisease.nf` that belongs in a subworkflow.
- **Reuse over duplication** — `RANK_VARIANTS`, `ANNOTATE_CSQ_PLI`, and `VCF_FILTER_BCFTOOLS_ENSEMBLVEP` are intentionally included multiple times under different aliases. Follow this pattern before creating a near-identical subworkflow.
- **nf-core modules first** — prefer a module from `modules/nf-core/` over writing a local one. Only add to `modules/local/` when no nf-core module exists or the use case is too pipeline-specific.

### Adding a new step

1. Define the corresponding input channel into your new process from the expected previous process channel.
2. Write the process block.
3. Define the output channel if needed.
4. Add any new parameters to `nextflow.config` with a default.
5. Add any new parameters to `nextflow_schema.json` with help text (via the `nf-core pipelines schema build` tool).
6. Add sanity checks and validation for all relevant parameters.
7. Perform local tests to validate that the new code works as expected.
8. If applicable, add a new test in the `tests` directory.
9. Update MultiQC config `assets/multiqc_config.yml` so relevant suffixes, file name clean up and module plots are in the appropriate order. If applicable, add a [MultiQC](https://multiqc.info/) module.
10. Add a description of the output files and if relevant any appropriate images from the MultiQC report to `docs/output.md`.

### Channel conventions

- **Skip lists**: use `parseSkipList()` and the `skip_tools` / `skip_subworkflows` params. Don't gate logic with raw string comparisons against params directly.
- **Conditional channels**: always initialize to `channel.empty()` before any `if` block that may or may not assign them. Never leave a channel potentially undefined.
- **Channel helpers**: use `channelFromPath`, `channelFromPathWithMeta`, and `channelFromSamplesheet` from `utils_nfcore_raredisease_pipeline` rather than rolling your own `channel.fromPath` calls.

### Params & analysis types

- `params` must only be accessed in the main unnamed workflow (`workflow` in `main.nf`). Subworkflows and named workflows receive all values as explicit `val_*` arguments. Never reference `params` directly inside a subworkflow.
- New tools that only apply to `wgs`, `wes`, or `mito` must be gated on `val_analysis_type`.
- Skippable tools must be added to the `--skip_tools` or `--skip_subworkflows` param and handled via `parseSkipList()`.
- Default params go in `nextflow.config`. Don't hardcode values that a user might reasonably want to change. Once added, run `nf-core pipelines schema build` to register them in `nextflow_schema.json`.

### Publishing

The pipeline uses Nextflow's `publish:` block and `output {}` API for file publishing. Each subworkflow exposes its outputs as named typed channel emits; the top-level `publish:` block in `main.nf` mixes them into destination-named entries.

- Emit every publishable output as its own named typed channel — one emit per file type, no `ch_publish` tuple wrapping and no grouped mix inside the subworkflow.
- In `main.nf`, mix all channels that share a destination into **one** `publish:` entry and **one** `output {}` entry. The mixing belongs at the routing layer, not inside the subworkflow.
- Channels consumed by downstream processes (e.g. MultiQC) and also published are emitted once; the caller wires the same channel to both consumers.

#### Emit naming convention

Use `<process_or_alias>_<emit_name>` (lowercase, underscored) inside the subworkflow's `emit:` block:

- Use the **alias name** as the prefix when a process is imported with `as` — the alias already encodes the distinction (e.g. `PICARD_COLLECTWGSMETRICS as PICARD_COLLECTWGSMETRICS_WG` → prefix `picard_collectwgsmetrics_wg`).
- Append the **module's emit name** verbatim.
- Drop obvious redundancy when the emit name exactly repeats a word already in the process/alias name (e.g. `sentieon_wgsmetrics_wg_wgs_metrics` → `sentieon_wgsmetrics_wg_metrics`). Do not rename to describe the file format — always use the emit name.
- For `VERIFYBAMID_VERIFYBAMID2`, drop the repetition: use prefix `verifybamid_`.

| Layer                     | Convention                                                          | Example                                                               |
| ------------------------- | ------------------------------------------------------------------- | --------------------------------------------------------------------- |
| Subworkflow `emit:`       | `<process_or_alias>_<emit_name>`                                    | `mosdepth_global_txt`                                                 |
| `raredisease.nf` variable | `ch_<subworkflow>_<semantic_suffix>`                                | `ch_qc_bam_mosdepth_global_txt`                                       |
| `NFCORE_RAREDISEASE` emit | `<subworkflow>_<emit_name>`                                         | `qc_bam_mosdepth_global_txt`                                          |
| `publish:` entry          | one entry per destination, mixing all channels for that destination | `qc_bam = NFCORE_RAREDISEASE.out.qc_bam_mosdepth_global_txt.mix(...)` |

The **semantic suffix** is the part of the emit name that describes what the data is, not which tool produced it. When the subworkflow emit name starts with a process/module name, drop that prefix in the `raredisease.nf` variable if the remainder is unambiguous within the subworkflow's outputs:

- `scatter_genome` emits `gatk4_splitintervals_split_intervals` → variable is `ch_scatter_genome_split_intervals` (drop `gatk4_splitintervals_`)
- `qc_bam` emits `mosdepth_global_txt` → variable stays `ch_qc_bam_mosdepth_global_txt` (`global_txt` alone would be ambiguous among the many txt outputs in that subworkflow)

When in doubt, keep enough of the process name to remain unambiguous.

> **Note:** Some subworkflows still use the legacy `ch_publish`/`subworkflow_results` pattern and are being migrated incrementally. Until a subworkflow is migrated, follow the existing pattern for that subworkflow so it continues to publish correctly via `subworkflow_results`.

### Configuration

- Process-level options go in `conf/modules/<subworkflow_name>.config`, not inline in the subworkflow `.nf` file.
- Only `ext.args`, `ext.args2`, and `ext.prefix` belong in module configs. Don't add business logic there.
- Conditional behavior (e.g. save as CRAM vs BAM) is handled in the subworkflow via `channel.empty()` gating — not via config-level flags.
- Process resource requirements (CPUs / memory / time) go in `conf/base.config` using `withLabel:` selectors so they can be shared across processes. Use `${task.cpus}` and `${task.memory}` in `script:` blocks to apply them dynamically.

### Writing tests

- Every subworkflow should have a test at `subworkflows/local/<name>/tests/main.nf.test`.
- Use `-stub` in the `when:` block only when real test data is difficult to generate. Prefer running with real data where it is reasonably available.
- Snapshot files (`*.nf.test.snap`) are committed alongside tests — update them when outputs change.
- Pipeline-level tests live in `tests/` and cover `default`, `test_bam`, and `test_singleton` profiles.
- Run `nf-test test <path>` for a single test, `nf-test test` for all.

### Style

- Sort `include` statements alphabetically by the name inside the braces. Right-pad each name with spaces so all closing `}` align to the same column (the longest name in the block sets the width):

  ```groovy
  include { ALIGN_BWA_BWAMEM2_BWAMEME                  } from '../align_bwa_bwamem2_bwameme'
  include { ALIGN_MT                                   } from '../align_MT'
  include { ALIGN_MT as ALIGN_MT_SHIFT                 } from '../align_MT'
  include { SAMTOOLS_VIEW as CONVERTTOCRAM_ALTFILTERED } from '../../../modules/nf-core/samtools/view/main'
  include { SAMTOOLS_VIEW as CONVERTTOCRAM_UNFILTERED  } from '../../../modules/nf-core/samtools/view/main'
  include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_EXCLUDE_ALT } from '../../../modules/nf-core/samtools/view/main'
  ```

- Both `take:` and `emit:` block entries require an inline type comment. Use `name // type: [mandatory|optional] description` for `take:` and `name = value // channel: [type description]` for `emit:`. Always include the comment — never leave an entry uncommented.

  ```groovy
  take:
      ch_vcf                // channel: [mandatory] [ val(meta), path(vcf) ]
      ch_reduced_penetrance // channel: [optional]  [ path(penetrance) ]
      val_aligner           // string:  [mandatory] aligner name (bwa/bwamem2/bwameme)
      process_with_sort     // Boolean

  emit:
      vcf = ch_vcf // channel: [ val(meta), path(vcf) ]
  ```

### Adding citations

When adding a new tool to the pipeline, update the following three locations:

#### 1. `CITATIONS.md`

Add an entry for the tool in alphabetical order under `## Pipeline tools`. If the tool has a publication, include a `>` citation block:

```markdown
- [ToolName](https://link-to-paper-or-repo)

  > Author A, Author B. Title. Journal. Year;vol(issue):pages. doi:...
```

If the tool has no publication, list only the link:

```markdown
- [ToolName](https://github.com/org/toolname)
```

#### 2. `subworkflows/local/utils_nfcore_raredisease_pipeline/main.nf`

Add citation text and bibliography entries inside `toolCitationText()` and `toolBibliographyText()`. Both functions are structured identically — group the tool's entry under the relevant category variable (e.g. `align_text`, `qc_bam_text`, `preprocessing_text`, `snv_annotation_text`). Mirror any conditional logic that gates the tool's execution (e.g. skip params, analysis type, or input content) so the citation only appears when the tool actually runs:

```groovy
// toolCitationText()
qc_bam_text = [
    ...,
    (condition) ? "ToolName (Author et al., Year)," : ""
]

// toolBibliographyText()
qc_bam_text = [
    ...,
    (condition) ? "<li>Author A, Author B. Title. Journal. Year. doi:...</li>" : ""
]
```

For tools that run only when the input samplesheet contains a particular file type, use a helper function rather than a param check — see `hasSpringInput()` as an example.

#### 3. `README.md`

Add the tool to the relevant numbered section in the **Pipeline summary**. If the tool belongs to a new category not yet represented, add a new numbered section in the appropriate position.
