# AREE methylation / WGBS / EM-seq workflow

STATUS: **unexecuted structural scaffold.** This directory contains real
DSL2 Nextflow wiring and real per-tool command lines, but no process here has
been run in this build — there is no compute budget or real/synthetic
bisulfite FASTQ data in this repository. Read this file before assuming
anything below has been validated end-to-end.

## What this workflow does

Implements Layer 2.B of the AREE design (`CLAUDE.md`, "DNA methylation / WGBS
/ EM-seq"): read QC, alignment, methylation calling, coverage filtering,
DML/DMR analysis, annotation to genomic context, and a standardized
genomic-region result table — or, when raw data are unavailable, harmonizing
an already-processed DMR table into the same standardized schema. Both modes
converge on the identical output schema and both terminate in an emitted
manifest and an HTML report, per `docs/design.md` Sec. 7.

## The two modes

Controlled by `params.mode`, default `processed_results_harmonization`.

### `raw_reanalysis`

For studies where raw bisulfite/EM-seq FASTQ are available. Pipeline:

```
FASTQC ─┐
        ├─> (QC only, does not block downstream)
TRIM_GALORE
   │
   ▼
BISMARK_GENOME_PREPARATION (once per reference genome)
   │
   ▼
BISMARK_ALIGN ─> BISMARK_DEDUPLICATE ─> BISMARK_METHYLATION_EXTRACTOR
   │
   ▼
COVERAGE_FILTER (per sample)
   │
   ▼
DMR_METHYLKIT (methylKit: methRead -> filterByCoverage -> unite ->
               calculateDiffMeth -> getMethylDiff, pooled across samples)
   │
   ▼
ANNOTATE_REGIONS (GenomicRanges/rtracklayer overlap against a GTF:
                  promoter / exon / intron / gene_body / intergenic)
   │
   ▼
STANDARDIZE_OUTPUT (reshape to the exact 9-column schema below)
   │
   ▼
EMIT_MANIFEST  ─>  RENDER_REPORT
```

Required inputs:

- `params.reads` — glob for paired-end FASTQ, e.g.
  `'data/raw/methylation/*_R{1,2}.fastq.gz'`
- `params.genome_fasta` — reference genome FASTA
- `params.annotation_gtf` — reference GTF (must have `gene`, `exon`, and
  ideally `transcript` feature rows; falls back to `gene` rows as transcript
  proxies if `transcript` rows are absent)
- `params.methylation.sample_sheet` — CSV with header `sample_id,treatment`
  (`treatment` is `0` = control/sham, `1` = treatment/exposed), used to build
  methylKit's `sample.id`/`treatment` vectors
- `params.methylation.study_id`, `params.methylation.comparison_id`

**This mode has never been run.** There is no synthetic bisulfite FASTQ
fixture in this repository to smoke-test it against (see
`config/demo.config` header). A future contributor adding a raw-data smoke
test should generate a small synthetic bisulfite-converted FASTQ pair (e.g.
via `Sherman` or a hand-rolled converter) plus a toy genome/GTF, per the
roadmap in `docs/roadmap.md`.

### `processed_results_harmonization`

For studies where only a processed DMR/DML table is publicly available (the
common case for older or third-party-reanalyzed methylation studies). Pipeline:

```
STANDARDIZE_OUTPUT (validate/reshape the supplied DMR table)
   │
   ▼
EMIT_MANIFEST (adds a warning: raw QC/coverage/DMR-calling params were
               not re-run or independently verified)
   │
   ▼
RENDER_REPORT
```

Required inputs:

- `params.methylation.processed_results` — path to a DMR-shaped TSV/CSV. At
  minimum needs `chrom`, `start`, `end` columns; `region_id`, `gene_id`,
  `annotation_context`, `meth_diff_percent` (or `meth_diff`), `qvalue` (or
  `q_value`), and `direction` are used when present and re-derived/left blank
  when absent (direction is inferred from the sign of the methylation
  difference if not supplied).
- `params.methylation.study_id`, `params.methylation.comparison_id`

`config/demo.config` wires this mode against the repository's demo DMR table
(`data/demo/methylation/GIGAS_PATH03_OsHV1_challenge_vs_sham_dmr_demo.tsv`),
which is clearly labeled simulated data.

## Standardized output schema

`STANDARDIZE_OUTPUT` emits `<study_id>_<comparison_id>_dmr_standardized.tsv`
with exactly these tab-separated columns, matching
`data/demo/methylation/GIGAS_PATH03_OsHV1_challenge_vs_sham_dmr_demo.tsv` and
what `src/harmonize/methylation.py` expects:

| Column | Notes |
|---|---|
| `region_id` | source-assigned or methylKit-assigned region identifier |
| `chrom` | chromosome/scaffold |
| `start` | region start coordinate |
| `end` | region end coordinate |
| `gene_id` | **optional/nullable.** Intergenic regions have no `gene_id` and are *kept*, never dropped (see `src/harmonize/methylation.py`, which routes gene-less rows to `feature_type = genomic_region` with `mapping_confidence = unresolved` rather than discarding them) |
| `annotation_context` | `promoter` \| `exon` \| `intron` \| `gene_body` \| `intergenic` |
| `meth_diff_percent` | signed methylation difference, treatment − control |
| `qvalue` | multiple-testing-adjusted p-value from methylKit (or as supplied in processed mode) |
| `direction` | `hyper` \| `hypo`, derived from the sign of `meth_diff_percent` when not directly supplied |

The only rows genuinely dropped by `STANDARDIZE_OUTPUT` are ones missing
`chrom`/`start`/`end` entirely (no genomic coordinates at all) — that count
is reported in the per-run QC table and the manifest, never silently
discarded.

## Feeding into harmonization

Both modes' final standardized TSV is designed to be passed directly to:

```bash
aree harmonize --study STUDY_ID --input path/to/STUDY_ID_COMPARISON_dmr_standardized.tsv
```

which calls `src/harmonize/methylation.py::harmonize_methylation()` to
produce evidence records in the shared cross-study schema (Layer 3).

## Provenance and manifest

`EMIT_MANIFEST` writes `<study_id>_<comparison_id>_manifest.json` containing:
sha256 checksums of the standardized table and upstream QC files, the full
explicit parameter set used, hardcoded software versions matching the pinned
container tags in `containers/README.md` (FastQC 0.12.1, Trim Galore 0.6.10,
Bismark 0.24.2, Bioconductor `RELEASE_3_18` methylKit/GenomicRanges),
workflow name/version, `date_generated` (the Nextflow run's start time),
QC metrics (computed from the standardized table where feasible, `null` with
a documented reason otherwise), and a `warnings` array. In
`processed_results_harmonization` mode this array always includes an
explicit `raw_data_not_available` warning, per `docs/design.md` Sec. 7 — this
is intentionally propagated to `quality_flags` at the evidence-record level
downstream, not hidden.

A separate `versions.yml` (nf-core-style `process: {tool: version}` mapping)
is emitted by every process, using the actual CLI's own `--version` output
where the tool is genuinely installed in its declared container (not
hardcoded, unlike the manifest's `software_versions` block, which is
necessarily hardcoded since no container was pulled in this build).

## Containers

Every process uses the exact image tags declared in `containers/README.md`:
`biocontainers/fastqc:0.12.1--hdfd78af_0`,
`biocontainers/trim-galore:0.6.10--hdfd78af_0`,
`biocontainers/bismark:0.24.2--hdfd78af_1`,
`bioconductor/bioconductor_docker:RELEASE_3_18` (methylKit DMR calling,
GenomicRanges/rtracklayer annotation), `python:3.11-slim` (standardize/
manifest steps, stdlib only), `ghcr.io/quarto-dev/quarto:1.5.57` (report).
None of these images have been pulled or run as part of this task — see
`containers/README.md` "What is real vs. aspirational here."

## Running

```bash
# processed_results_harmonization mode against the shipped demo table
cd workflows/methylation
nextflow run main.nf -profile docker,demo

# raw_reanalysis mode (not runnable today — no FASTQ fixtures exist)
nextflow run main.nf -profile docker \
    --mode raw_reanalysis \
    --reads 'data/raw/methylation/*_R{1,2}.fastq.gz' \
    --genome_fasta data/reference/cgigas_genome.fa \
    --annotation_gtf data/reference/cgigas_annotation.gtf \
    --methylation.sample_sheet data/raw/methylation/sample_sheet.csv \
    --methylation.study_id STUDY_ID \
    --methylation.comparison_id COMPARISON_ID
```

**Neither of the above commands has actually been executed in this build.**
The `docker,demo` invocation is expected to be runnable on a machine with
Docker and Nextflow installed, since `processed_results_harmonization` mode
only depends on `python:3.11-slim` and `ghcr.io/quarto-dev/quarto:1.5.57`,
both lightweight, real, pullable images — but that has not been verified
here.

## What has and has not been run (honesty statement)

- **Not run:** every process in this workflow, in either mode.
- **Not run:** no container has been pulled or tested against these exact
  command lines.
- **Not verified:** the R scripts in `dmr_methylkit.nf` and
  `annotate_regions.nf` are syntactically correct, real methylKit/
  GenomicRanges idiom, written to the best of the author's knowledge of those
  APIs — they have not been executed against real or synthetic sequencing
  data, so subtle argument/behavior mismatches with a specific package
  version are possible and should be checked before first real use.
- **Real:** the CLI invocations for FastQC, Trim Galore, and all three
  Bismark steps are genuine, correctly-flagged commands for paired-end
  bisulfite data.
- **Real:** `STANDARDIZE_OUTPUT` and `EMIT_MANIFEST` use only the Python
  standard library and would execute correctly today if run directly with
  `python3` against the demo TSV — they were not, however, invoked as part of
  this task, and have not been run inside Nextflow.
- See also the repository-wide implementation-status table (`docs/roadmap.md`
  or the top-level status table referenced from `README.md`) for how this
  workflow's status compares to the other three assay workflows.
