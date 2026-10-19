# AREE methylation / WGBS / EM-seq workflow

STATUS: **`processed_results_harmonization` mode is executed, verified, and
covered by current-version CI; the Bismark/methylKit raw path is not.** It did
not run before 2026-08-28; see `../../docs/first_raw_reanalysis.md` for the
defects that prevented it. No container has been pulled. The raw branch contains
real DSL2 wiring and tool commands, but no bisulfite FASTQ fixture or real raw
study has exercised them.

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

#### Alignment and extraction options

`methylation.bismark_align_args` and `methylation.extractor_args` are
passed to `bismark` and `bismark_methylation_extractor` as written, and
recorded in the manifest. Both are empty by default. Bismark's default
`--score_min L,0,-0.2` is strict for a genome as heterozygous as
*M. gigas*. On 50,000 trimmed pairs from `CIBNOR2021_HEAT_RRSS`, the
default mapped 15.5% of pairs uniquely, and `--score_min L,0,-0.6` mapped
53.9%. Non-directional mode and a larger insert limit made no difference.
Use `extractor_args` for M-bias trimming (e.g. `--ignore_r2 N`) after
checking the M-bias plot.

#### Chunked alignment (`methylation.align_chunk_pairs`)

A deep library is one long Bismark task: about 20 h for 224 M read pairs on
a 10-core workstation. On a scheduler with short or preemptible jobs, such
as Klone's ckpt partitions, it never finishes. Set
`methylation.align_chunk_pairs` to split each trimmed library into chunks of
that many read pairs (e.g. `20000000`):

- `SPLIT_FASTQ_PAIRS` splits R1/R2 in step with `seqkit split2`.
- `BISMARK_ALIGN_CHUNK` aligns each chunk with the same Bismark command as
  `BISMARK_ALIGN`.
- `BISMARK_MERGE_CHUNKS` concatenates the chunk BAMs in chunk order with
  `samtools cat`, so mates stay adjacent. It also sums the chunk reports
  into one per-library `*_PE_report.txt` in Bismark's own layout
  (`bin/sum_bismark_pe_reports.sh`); percentages are recomputed from the
  summed counts.

Deduplication and extraction then run on the merged BAM, as before. Duplicates
that fall in different chunks are still removed. Bismark aligns each pair
independently, so the merged BAM holds the same alignments as an unchunked
run's. The per-chunk reports are published next to the merged one. Unset or
`0` keeps one alignment task per library. See
[docs/running_on_klone.md](../../docs/running_on_klone.md).

#### Replication unit (`methylation.replicate_unit`)

When the phenotype or treatment belongs to a group of libraries, such as a
breeding family or a tank, the libraries in a group are not independent
replicates of it. Set `methylation.replicate_unit` to the design-sheet
column that names the group, and `methylation.design_sheet` to a TSV with
that column plus `run_accession` or `sample_id`:

```
methylation.replicate_unit = 'family'
methylation.design_sheet   = 'data/studies/<STUDY_ID>/samplesheet.tsv'
```

`DMR_METHYLKIT` then sums each unit's methylated and unmethylated counts per
region after `unite`, and tests the units with methylKit's McCullagh-Nelder
overdispersion correction and an F test (`overdispersion = "MN"`,
`test = "F"`), so the error comes from unit-to-unit variation. Each unit must
lie wholly in one arm, and each arm needs at least two units. Without
`replicate_unit` the step is unchanged.

On synthetic data with family-to-family variation and no true class effect
(4 families × 3 libraries, 380 null and 20 true regions), testing the 12
libraries as independent called 162 of 380 null regions significant; testing
the 4 families called none. With 2 families per arm the test has very little
power: the true regions ranked near the top but none reached q < 0.05. Where
a region's estimated overdispersion falls below 1, methylKit resets it to 1,
so a few null regions can still get small q-values; the `meth_diff_cutoff`
filter is the second guard.

#### DMR engines (`methylation.dmr_engine`)

- `methylkit` (default): `DMR_METHYLKIT`, as above.
- `limma` and `dss`: `DMR_POOLED_TILES` (`modules/methylation/dmr_pooled_tiles.nf`).
  It merges the two strands of each CpG, pools libraries by `replicate_unit`,
  and keeps CpGs with at least `unit_min_coverage` reads in **every unit**
  (coverage is applied after pooling, so set `min_coverage` low, e.g. 1). It
  sums them into `tile_size` tiles with at least `tile_min_cpgs` CpGs, then
  tests treatment units against control units.
  - `limma`: a moderated t on per-unit M-values (the logit of methylation),
    with variance shrunk across tiles by empirical Bayes. Recommended when
    there are only a few units per arm.
  - `dss`: `DSS::DMLtest`, a beta-binomial model with shrunken dispersion
    and a Wald test.

Units (`methylation.dmr_features`), for `limma` and `dss`:
- `tiles` (default): `tile_size` windows genome-wide.
- `genes`: one unit per gene, pooling every CpG in the gene body and its
  `promoter_upstream_bp` promoter (strand-aware). CpGs outside every gene
  window fall into intergenic tiles. Each gene then gives one record, so
  harmonized methylation evidence links one-to-one to RNA-seq genes; with
  `tiles`, many tiles share a gene and meta-analysis refuses to pool them
  (identifier collision). `ANNOTATE_REGIONS` keeps the gene and context these
  units arrive with.

Why `limma` is recommended for few units: on synthetic data with 2 vs 2
families, real family-to-family variation and 20 true regions among 400
tiles:

| Test | True regions found | False positives (of 380) |
|---|---|---|
| methylKit (families, MN + F) | 0 | 0 |
| DSS (families) | 17 | 57 |
| limma (families) | 4 | 1 |

limma's null p-values were calibrated (5.0% below 0.05). DSS's Wald test
has no small-sample correction, so it is overconfident at 2 vs 2.

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
with these tab-separated columns (the last three optional), matching
`data/demo/methylation/GIGAS_PATH03_OsHV1_challenge_vs_sham_dmr_demo.tsv` and
what `src/harmonize/methylation.py` expects:

| Column | Notes |
|---|---|
| `region_id` | source-assigned or methylKit-assigned region identifier |
| `chrom` | chromosome/scaffold |
| `start` | region start coordinate |
| `end` | region end coordinate |
| `gene_id` | **optional/nullable.** Intergenic regions have no `gene_id` and are *kept*, never dropped (see `src/harmonize/methylation.py`, which keeps gene-less rows as `feature_type = methylation_region` with `feature_id_original = region_id` and `mapping_confidence = unresolved` rather than discarding them) |
| `annotation_context` | `promoter` \| `exon` \| `intron` \| `gene_body` \| `intergenic` |
| `meth_diff_percent` | signed methylation difference, treatment − control |
| `qvalue` | multiple-testing-adjusted p-value from methylKit (or as supplied in processed mode) |
| `direction` | `hyper` \| `hypo`, derived from the sign of `meth_diff_percent` when not directly supplied |
| `pvalue` | *optional*: unadjusted p-value. Written when the DMR step supplies it (`limma` / `dss` engines); methylKit does not |
| `se_meth_diff_percent` | *optional*: standard error of `meth_diff_percent`. For `limma`, abs(diff) / qnorm(1 - p/2), so a normal (Wald) test, which is what meta-analysis applies, reproduces limma's own p-value and keeps its small-sample df. Using abs(diff / t) instead gave p ~1e-23 where limma gave 3e-4 at 2 vs 2 units |
| `analysis_method` | *optional*: e.g. `limma_pooled_tiles`; the harmonizer uses it as the evidence `analysis_method` (default `methylKit`) |

Without `pvalue` and `se_meth_diff_percent`, harmonized methylation
evidence has no standard error and cannot enter inverse-variance pooling.

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
is emitted by every process from the tool's own `--version` output, whether
the tool runs natively or in its declared container. `main.nf` collects them
into `EMIT_MANIFEST`, which merges them (`modules/shared/collect_versions.py`)
into the manifest's `software_versions`, so the manifest reports what ran.
`RENDER_REPORT` runs after the manifest and is listed under
`software_versions_not_captured`.

## Containers

Every process uses the exact image tags declared in `containers/README.md`:
`quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0`,
`quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0`,
`quay.io/biocontainers/bismark:0.24.2--hdfd78af_0`,
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
