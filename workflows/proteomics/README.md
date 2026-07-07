# AREE proteomics workflow

STATUS: **unexecuted structural scaffold.** This directory contains real
DSL2 Nextflow wiring and real per-tool logic, but no process here has been
run in this build — there is no compute budget or real mass-spec data in
this repository. Read this file before assuming anything below has been
validated end-to-end.

## What this workflow does

Implements Layer 2.C of the AREE design (`CLAUDE.md`, "Proteomics"): input
harmonization for protein abundance tables, normalization, a missingness
report, differential abundance, protein-to-gene identifier translation, and a
standardized protein-level evidence output — or, when only a processed
differential-abundance table is available, harmonizing that table directly
into the same standardized schema. Both modes converge on the identical
output schema and both terminate in an emitted manifest and an HTML report,
per `docs/design.md` Sec. 7.

## The two modes

Controlled by `params.mode`, default `processed_results_harmonization`.

### `raw_reanalysis`

For studies that deposited a wide, already-quantified protein abundance
matrix (protein_accession x sample). **"Raw" here means the pre-normalization
abundance matrix, not raw mass-spec spectra (.raw/.mzML)** — full spectral
search (MaxQuant/FragPipe-level) is out of scope for this scaffold; see
`docs/roadmap.md`. Pipeline:

```
HARMONIZE_INPUT (validate matrix shape, sample sheet, protein_accession column)
   │
   ▼
NORMALIZE (median/quantile normalization across samples)
   │
   ▼
MISSINGNESS_REPORT (per-protein and per-sample missingness before any
                     imputation; reported, never silently imputed away)
   │
   ▼
DIFFERENTIAL_ABUNDANCE (limma, moderated t-statistics)
   │
   ▼
ID_TRANSLATION (protein_accession -> gene identifier via
                data/mappings/gene_id_crosswalk.tsv)
   │
   ▼
STANDARDIZE_OUTPUT (reshape to the exact schema below)
   │
   ▼
EMIT_MANIFEST  ─>  RENDER_REPORT
```

Required inputs:

- `params.proteomics.raw_abundance_matrix` — wide TSV/CSV,
  `protein_accession` column + one column per sample
- `params.proteomics.sample_sheet` — TSV/CSV with `sample_id`, `group` (+
  optional metadata columns)
- `params.proteomics.id_mapping_table` — defaults to
  `data/mappings/gene_id_crosswalk.tsv`
- `params.proteomics.study_id`, `params.proteomics.comparison_id`

**This mode has never been run.** There is no synthetic raw abundance-matrix
fixture in this repository to smoke-test it against — `data/demo/proteomics/`
only ships an already-differential result table (see `config/demo.config`).
A future contributor adding a raw-data smoke test should generate a small
synthetic wide abundance matrix, per the roadmap in `docs/roadmap.md`.

### `processed_results_harmonization`

For studies that only released a processed protein-level differential
abundance table (the common case — see `registry/studies/GIGAS_SAL04.yaml`,
which is exactly this scenario). Pipeline:

```
STANDARDIZE_OUTPUT (validate/reshape the supplied table)
   │
   ▼
EMIT_MANIFEST (adds a warning: raw QC/normalization/missingness-handling
               were not re-run or independently verified)
   │
   ▼
RENDER_REPORT
```

Required inputs:

- `params.proteomics.processed_results` — path to a protein-level TSV/CSV
  matching the schema below
- `params.proteomics.study_id`, `params.proteomics.comparison_id`

`config/demo.config` wires this mode against the repository's demo protein
abundance table
(`data/demo/proteomics/GIGAS_SAL04_low_salinity_vs_control_protein_abundance_demo.tsv`),
which is clearly labeled simulated data. When no independently-computed
missingness report exists for this mode, the workflow passes the sentinel
file `assets/NO_FILE` through instead of a real path — `standardize_output.nf`
and `emit_manifest.nf` check for that literal filename and record
"not supplied" rather than treating it as a real, empty result.

## Standardized output schema

`STANDARDIZE_OUTPUT` emits `<study_id>_<comparison_id>_protein_standardized.tsv`
with exactly these tab-separated columns, matching
`data/demo/proteomics/GIGAS_SAL04_low_salinity_vs_control_protein_abundance_demo.tsv`
and what `src/harmonize/proteomics.py` expects:

| Column | Notes |
|---|---|
| `protein_accession` | UniProt accession (or source-native protein ID) |
| `gene_symbol` | as reported by the source study; not used for standardization, informational only — `src/harmonize/proteomics.py` resolves identity via `protein_accession` against the crosswalk |
| `log2FC` | signed log2 fold-change, treatment vs control |
| `pvalue` | raw p-value from the differential abundance step |
| `padj` | multiple-testing-adjusted p-value |
| `missingness_percent` | percent of samples missing a value for this protein prior to any imputation |

## Feeding into harmonization

Both modes' final standardized TSV is designed to be passed directly to:

```bash
aree harmonize --study STUDY_ID --input path/to/STUDY_ID_COMPARISON_protein_standardized.tsv
```

which calls `src/harmonize/proteomics.py::harmonize_proteomics()` to produce
evidence records in the shared cross-study schema (Layer 3), including a
`processed_only` quality flag whenever `missingness_percent >= 10` or the
study's own registration already declares processed-only data availability.

## Provenance and manifest

`EMIT_MANIFEST` writes `<study_id>_<comparison_id>_manifest.json` containing:
sha256 checksums of the standardized table and upstream QC files, the full
explicit parameter set used, hardcoded software versions matching the pinned
container tags in `containers/README.md`, workflow name/version,
`date_generated` (the Nextflow run's start time), QC metrics (missingness
summary where available, `null` with a documented reason otherwise), and a
`warnings` array. In `processed_results_harmonization` mode this array always
includes an explicit `raw_data_not_available` warning, propagated to
`quality_flags` at the evidence-record level downstream, not hidden.

## Containers

Every process uses the exact image tags declared in `containers/README.md`:
`python:3.11-slim` (input harmonization, normalization, missingness report,
ID translation, standardize/manifest steps), `bioconductor/bioconductor_docker:RELEASE_3_18`
(limma differential abundance), `ghcr.io/quarto-dev/quarto:1.5.57` (report).
None of these images have been pulled or run as part of this task — see
`containers/README.md` "What is real vs. aspirational here."

## Running

```bash
# processed_results_harmonization mode against the shipped demo table
cd workflows/proteomics
nextflow run main.nf -profile docker,demo

# raw_reanalysis mode (not runnable today — no abundance-matrix fixture exists)
nextflow run main.nf -profile docker \
    --mode raw_reanalysis \
    --proteomics.raw_abundance_matrix data/raw/proteomics/abundance_matrix.tsv \
    --proteomics.sample_sheet data/raw/proteomics/sample_sheet.csv \
    --proteomics.study_id STUDY_ID \
    --proteomics.comparison_id COMPARISON_ID
```

**Neither of the above commands has actually been executed in this build.**

## What has and has not been run (honesty statement)

- **Not run:** every process in this workflow, in either mode.
- **Not run:** no container has been pulled or tested against these exact
  command lines.
- **Not verified:** the limma R script in `differential_abundance.nf` is
  syntactically correct, real limma idiom, written to the best of the
  author's knowledge of that API — it has not been executed against real or
  synthetic proteomics data, so subtle argument/behavior mismatches with a
  specific package version are possible and should be checked before first
  real use.
- **Real:** `HARMONIZE_INPUT`, `NORMALIZE`, `MISSINGNESS_REPORT`,
  `ID_TRANSLATION`, `STANDARDIZE_OUTPUT`, and `EMIT_MANIFEST` use only the
  Python standard library (plus pandas) and would execute correctly today if
  run directly against the demo TSV — they were not, however, invoked as
  part of this task, and have not been run inside Nextflow.
- See also the repository-wide implementation-status table (`docs/roadmap.md`
  or the top-level status table referenced from `README.md`) for how this
  workflow's status compares to the other three assay workflows.
