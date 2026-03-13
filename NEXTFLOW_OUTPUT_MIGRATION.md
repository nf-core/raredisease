---

## Centralized Publishing Pattern (nf-core migration approach)

Instead of using `topic` channels or `publishDir`, we use a **single `ch_publish` channel
per subworkflow** that emits `[destination, value]` tuples. Publishing is handled centrally
in `main.nf`.

### Pattern Overview

**Inside each subworkflow**, map every published output into a `[destination, value]` tuple
where `destination` mirrors the original `publishDir` path. Mix all of them into a single
`ch_publish` and emit it:

```nextflow
// Inside subworkflow
ch_publish = PROCESS_A.out.bam
    .map { meta, bam -> ['alignment/', [meta, bam]] }
    .mix(PROCESS_B.out.metrics.map { meta, f -> ['qc/', [meta, f]] })
    .mix(PROCESS_C.out.vcf.map { meta, vcf -> ['vcf/', [meta, vcf]] })

emit:
ch_publish
```

Rules:

- There must be only **one `ch_publish` emit per subworkflow**
- The destination string should mirror the original `publishDir` path
- Preserve the original channel structure as the second element of the tuple
- Include outputs from processes whose channels are not consumed downstream —
  these were previously captured by `publishDir` at the filesystem level

**In `main.nf`**, collect `ch_publish` from all subworkflows, mix them, and publish
centrally:

```nextflow
workflow {
    main:
    SUBWORKFLOW_A(...)
    SUBWORKFLOW_B(...)

    publish:
    all_outputs = SUBWORKFLOW_A.out.ch_publish
        .mix(SUBWORKFLOW_B.out.ch_publish)
}

output {
    all_outputs {
        path { destination, value -> destination }
    }
}
```

### Handling varying channel structures

Channels with different numbers of elements (e.g. `[meta, bam]` vs `[meta, bam, bai]`)
are handled safely by always wrapping the entire inner tuple as a single second element.
Nextflow recursively scans the value for file objects automatically.

// CORRECT
.map { meta, bam, bai -> ['alignment/', [meta, bam, bai]] }

// WRONG — breaks the [destination, value] contract
.map { meta, bam, bai -> ['alignment/', meta, bam, bai] }
```

**Your command** — add one line to enforce the wrapping rule:
```
...

For each process output that should be published, map it into a tuple of
[destination_path, channel_value] where destination_path mirrors the directory
that publishDir was publishing to. Always wrap the entire channel value as a
single second element regardless of how many elements it contains, e.g.
{ meta, bam, bai -> ['alignment/', [meta, bam, bai]] }.
...

### Grouping channels by destination before mapping

Instead of mapping each channel individually, group channels sharing the same destination
using `mix` first, then apply a single `map` per destination group. This reduces repetition
significantly.

```nextflow
// PREFERRED — one map per destination group
ch_qc_bam = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
    .mix(PICARD_COLLECTMULTIPLEMETRICS.out.pdf)
    .mix(TIDDIT_COV.out.wig)
    .mix(MOSDEPTH.out.global_txt)
    .mix(ch_qualimap)
    .map { meta, value -> ['qc_bam/', [meta, value]] }

ch_ngsbits = ch_ngsbits_samplegender
    .map { meta, tsv -> ['ngsbits_samplegender/', [meta, tsv]] }

ch_publish = ch_qc_bam.mix(ch_ngsbits)

// AVOID — one map per channel, very repetitive
ch_publish = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
    .map { meta, metrics -> ['qc_bam/', [meta, metrics]] }
    .mix(PICARD_COLLECTMULTIPLEMETRICS.out.pdf
        .map { meta, pdf -> ['qc_bam/', [meta, pdf]] })
    // ...
```

Rules:

- Create one intermediate channel per unique destination directory
- Mix all channels going to the same destination first, then apply a single `map`
- Combine all destination groups into `ch_publish` with a final `mix`

### Nested subworkflows — bubbling up ch_publish

When a subworkflow calls other subworkflows, always mix the inner `ch_publish` into the
outer `ch_publish`. This propagates publishing up the call chain to `main.nf`.

```nextflow
// Inner subworkflow emits its own ch_publish
workflow ALIGN_BWA {
    main:
    BWA_MEM(...)
    SAMTOOLS_SORT(...)

    ch_publish = BWA_MEM.out.bam
        .mix(SAMTOOLS_SORT.out.bam)
        .map { meta, value -> ['alignment/', [meta, value]] }

    emit:
    ch_publish
}

// Outer subworkflow mixes in ch_publish from inner subworkflow
workflow ALIGN {
    main:
    ALIGN_BWA(...)
    FASTP(...)

    ch_publish = FASTP.out.reads
        .map { meta, value -> ['fastp/', [meta, value]] }
        .mix(ALIGN_BWA.out.ch_publish)

    emit:
    ch_publish
}
```

Rules:

- Every subworkflow that has publishable outputs must emit `ch_publish`
- If a subworkflow calls inner subworkflows, always mix their `ch_publish` into the
  outer `ch_publish` — never discard it
- In `main.nf`, mix `ch_publish` from all top-level subworkflows into one channel
  before the `publish:` block

### Why this approach

- Subworkflow signatures stay clean — one `ch_publish` emit regardless of how many
  output directories exist internally
- Publishing logic is fully centralized in `main.nf`
- No `topic` channel footguns (pipeline hanging forever if a process consumes its own topic)
- Easy to audit all published paths by reading `main.nf` alone
