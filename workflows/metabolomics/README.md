# AREE metabolomics workflow

STATUS: **both modes execute.** `processed_results_harmonization` is verified
against the bundled demo. `raw_reanalysis` runs natively (`-profile local`) on
the simulated targeted-panel fixture in `tests/fixtures/metabolomics_raw/`
(checked in CI by `scripts/check_metabolomics_raw_smoke.py`) and on the real
targeted panel of `DUPERRET2025_POMS_STARVE`. No container has been pulled.
See `../../docs/first_raw_reanalysis.md` for the defects that kept the
processed path from running before 2026-08-28.

This implements Layer 2.D ("Metabolomics") of the AREE reanalysis workflows
described in `CLAUDE.md`: metabolite feature table intake, annotation
confidence tracking, normalization and QC, differential abundance, pathway or
metabolite-class mapping, and standardized feature-level evidence output.

## Important simplification: what "raw" means here

Like the proteomics workflow, **"raw" in this workflow does not mean raw mass
spectrometry spectra (mzML/RAW files).** Reprocessing spectral data from
scratch — peak picking, chromatographic alignment, retention-time
correction — is a large, tool-specific undertaking (XCMS, MZmine, MS-DIAL,
etc.) that is explicitly **out of scope** for AREE and is not implemented
anywhere in this repository.

Instead, `raw_reanalysis` mode in this workflow takes as its raw input a
**raw/unfiltered metabolite feature intensity table** — the kind of
features-by-samples matrix that XCMS (`featureTable`) or MZmine (aligned
feature list export) produce after peak picking and alignment, before any
AREE-side filtering, normalization, or statistics have been applied. A lab
would run their preferred spectral-processing tool upstream and feed the
resulting feature table into `--metabolomics.raw_feature_table`.

This is a documented simplification, not an oversight — see
`docs/design.md` section 9 (Explicit assumptions) and the equivalent note in
`workflows/proteomics/README.md` for the same pattern applied to proteomics
raw abundance tables.

## Two modes

Set via `params.mode` (`--mode` on the CLI), matching `docs/design.md`
section 7:

### `raw_reanalysis`

```
params.metabolomics.raw_feature_table (features x samples intensities)
params.metabolomics.sample_sheet      (sample_id, condition[, qc_pool])
        │
        ▼
FEATURE_TABLE_INTAKE     — structural validation, numeric coercion
        │
        ▼
ANNOTATION_CONFIDENCE    — assign/validate MSI-style 1-4 confidence level
        │
        ▼
NORMALIZE_QC             — TIC normalization + log2; pooled-QC-sample CV
        │
        ▼
DIFFERENTIAL_ABUNDANCE   — limma moderated t-test (treatment vs control)
        │
        ▼
PATHWAY_MAPPING          — putative_metabolite_name -> metabolite_class/pathway
        │
        ▼
STANDARDIZE_OUTPUT       — reshape/validate to the exact 7-column schema
        │
        ├──▶ EMIT_MANIFEST   — checksums, params, versions, QC, warnings
        └──▶ RENDER_REPORT   — Quarto HTML report
```

### `processed_results_harmonization` (default)

```
params.metabolomics.processed_results (already-shaped feature table)
        │
        ▼
STANDARDIZE_OUTPUT       — VALIDATION GATE: reshape/validate to schema
        │
        ├──▶ EMIT_MANIFEST   — flags that raw QC/normalization were NOT
        │                      independently verified by this workflow
        └──▶ RENDER_REPORT   — Quarto HTML report
```

Both modes converge on the same `STANDARDIZE_OUTPUT` schema, so
`src/harmonize/metabolomics.py` never needs to know which mode produced a
given file — this mirrors the rnaseq/proteomics workflows and
`docs/design.md` section 7 exactly.

## Expected inputs

| Mode | Param | Description |
|---|---|---|
| raw_reanalysis | `--metabolomics.raw_feature_table` | Raw/unfiltered feature-by-sample intensity table (TSV/CSV). Must have a `feature_id` column plus one intensity column per `sample_id` in the sample sheet. May optionally already carry `putative_metabolite_name` / `annotation_confidence_level`. |
| raw_reanalysis | `--metabolomics.sample_sheet` | TSV/CSV with `sample_id`, `condition` (`control`/`treatment`), optional `qc_pool` (true/false) for pooled QC injections. |
| raw_reanalysis | `--metabolomics.metabolite_annotation_map` (optional) | TSV: `feature_id`, `putative_metabolite_name`, `annotation_confidence_level`. Used to fill in confidence levels not already present in the feature table. |
| raw_reanalysis | `--metabolomics.metabolite_class_map` (optional) | TSV: `putative_metabolite_name`, `metabolite_class`[, `pathway`]. See "Metabolite class/pathway crosswalk" below. |
| raw_reanalysis | `--metabolomics.qc_cv_threshold` | Coefficient-of-variation threshold for the pooled-QC-sample QC metric (default `0.30`). |
| raw_reanalysis | `--metabolomics.normalization` | `tic` (default), `median` or `internal_standard`; see "Targeted panels" below. |
| raw_reanalysis | `--metabolomics.detection_limit` | Values at or below this become missing (below detection), not intensities. Default `null`. |
| raw_reanalysis | `--metabolomics.replicate_unit` | Sample-sheet column (e.g. `tank`) whose units are the replicates; log2 values are averaged per unit before limma. Default `null`. |
| raw_reanalysis | `--metabolomics.min_detected_per_arm` | Test a feature only if detected in at least this many replicates per arm (default `2`). |
| raw_reanalysis | `--metabolomics.acquisition` | `untargeted` (default) or `targeted`; written to each row's `analysis_method`. |
| processed_results_harmonization | `--metabolomics.processed_results` | Already-shaped differential feature table (see target schema below), or a close variant `STANDARDIZE_OUTPUT` can coerce (tolerated header synonyms are listed in `modules/metabolomics/standardize_output.nf`). |
| both | `--metabolomics.study_id`, `--metabolomics.comparison_id` | Identifiers threaded into every output filename and the provenance manifest. |

## Targeted panels

A targeted MRM/SIM panel differs from an untargeted feature table in ways the
default TIC path ignores, so `raw_reanalysis` has options for it:

- **Platforms.** An optional `platform` column in the feature table (e.g. one
  value per LC or GC method) makes every normalization run within a platform,
  so panels on different instruments are never scaled against each other.
- **Internal standards.** Rows with `is_internal_standard` true are used by
  `--metabolomics.normalization internal_standard` (each sample's mean log2
  standard area, per platform, is subtracted) and are never tested.
  `median` centres each sample on its median detected metabolite instead.
- **Detection limit.** Vendor exports often write a placeholder such as
  `0.0001` below the limit of detection. `--metabolomics.detection_limit`
  turns those into missing values. Features detected in fewer than
  `--metabolomics.min_detected_per_arm` replicates of either arm are not
  tested, and `<study>_<comparison>_missingness.tsv` (in
  `differential_abundance/`) lists every feature with its detection counts
  and status, so nothing disappears silently.
- **Replicate units.** `--metabolomics.replicate_unit` averages samples that
  share a unit (e.g. two families pooled from one tank) before the fit, the
  metabolomics counterpart of the RNA-seq and methylation option. A unit may
  not sit in both arms.

Raw mode also writes `lfcSE` (limma's posterior standard error) and
`analysis_method`, which `aree harmonize` carries into the evidence table.

## Metabolite class/pathway crosswalk

`data/mappings/` currently contains only `gene_id_crosswalk.tsv` (a
**gene**-identifier crosswalk for NCBI/Ensembl/UniProt/locus/symbol/
orthogroup, used by RNA-seq and proteomics harmonization) and
`ambiguous_symbol_map.yaml`. Neither applies to metabolite names or classes,
and no metabolite-class/pathway crosswalk currently ships in this
repository. `PATHWAY_MAPPING` therefore requires an explicit
`--metabolomics.metabolite_class_map` TSV; if omitted, every feature is
conservatively assigned `metabolite_class = "unknown"` rather than guessed.
Adding a small illustrative `data/mappings/metabolite_class_map.tsv` (e.g.
mapping the metabolite names already used in
`data/demo/metabolomics/*_demo.tsv` to classes like `organic_acid`,
`amino_acid`, `nucleotide`, `osmolyte`, `carbohydrate`) is listed as
follow-up work, not implemented in this change (this task was scoped to
`workflows/metabolomics/` and `modules/metabolomics/` only).

## Target output schema

`STANDARDIZE_OUTPUT` emits `<study_id>_<comparison_id>_metabolite_features_standardized.tsv`,
tab-separated, with exactly these columns (matching
`data/demo/metabolomics/GIGAS_GROW06_hypoxia_vs_control_metabolite_features_demo.tsv`):

```
feature_id  putative_metabolite_name  annotation_confidence_level  metabolite_class  log2FC  pvalue  padj
```

plus, when the input has them (raw mode always does), `lfcSE` and
`analysis_method`.

`annotation_confidence_level` uses a simplified Metabolomics Standards
Initiative (MSI)-style 1-4 scale:

| Level | Meaning |
|---|---|
| 1 | Confirmed by authentic chemical standard (RT + MS/MS match) |
| 2 | Putatively annotated (spectral library match, no standard run) |
| 3 | Putative class / substructure only |
| 4 | Unknown / unannotated feature |

This is the exact contract expected by `src/harmonize/metabolomics.py`,
which maps `annotation_confidence_level` to the evidence-record
`mapping_confidence` field (`1→exact`, `2→inferred`, `3`/`4→unresolved`) —
see that file's module docstring for the rationale (only level 1-2 features
get a standardized identifier at all).

## Mapping workflow output to `aree harmonize`

Once `STANDARDIZE_OUTPUT` has produced the standardized TSV (in either
mode), it is a direct input to the Python harmonization CLI:

```bash
aree harmonize \
  --study GIGAS_GROW06 \
  --input results/demo/metabolomics/standardized/GIGAS_GROW06_hypoxia_vs_control_metabolite_features_standardized.tsv
```

which calls `src/harmonize/metabolomics.py::harmonize_metabolomics()` to
produce evidence records in the shared cross-study evidence schema
(`schemas/evidence.schema.json`).

## What has and has not been run

**Not run in this build:**
- `nextflow run workflows/metabolomics/main.nf` (any mode, any profile) —
  no Nextflow, Docker/Apptainer runtime, or compute budget available here.
- None of the `container` images listed below have been pulled or verified
  against these exact process command lines.
- No raw/unfiltered feature intensity table (real or synthetic) exists in
  this repository, so `raw_reanalysis` mode has never been exercised
  end-to-end, only reviewed for syntactic/structural correctness.

**What is real:**
- Every Python `script:` block is real pandas logic that reads/writes actual
  files with actual column-name validation, not placeholder echo statements.
- The R differential-abundance script is real, syntactically valid limma
  usage (`lmFit` → `eBayes` → `topTable`).
- The Quarto report template (`assets/report_template.qmd`) is a real
  parameterized R Markdown/Quarto document using `ggplot2`/`dplyr`/`knitr`.
- `config/demo.config` already wires `params.metabolomics.processed_results`
  to `data/demo/metabolomics/GIGAS_GROW06_hypoxia_vs_control_metabolite_features_demo.tsv`,
  so `-profile demo` (or `-c ../../config/demo.config`) demonstrates correct
  wiring of `processed_results_harmonization` mode, if actually run.

## Running (once a real environment is available)

```bash
# processed_results_harmonization mode against the repo's demo table
nextflow run workflows/metabolomics/main.nf \
  -profile docker,demo

# raw_reanalysis mode against your own feature table
nextflow run workflows/metabolomics/main.nf \
  -profile docker \
  --mode raw_reanalysis \
  --metabolomics.raw_feature_table path/to/features.tsv \
  --metabolomics.sample_sheet path/to/samples.tsv \
  --metabolomics.metabolite_class_map path/to/class_map.tsv \
  --metabolomics.study_id STUDY_ID \
  --metabolomics.comparison_id treatment_vs_control
```

## Container images used

See `containers/README.md` for the full cross-workflow strategy and honesty
notes. This workflow declares:

| Step | Image |
|---|---|
| Feature table intake / annotation confidence / normalization+QC / pathway mapping / standardize output / manifest (Python) | `python:3.11-slim` |
| Differential abundance (R, limma) | `bioconductor/bioconductor_docker:RELEASE_3_18` |
| Report render (Quarto) | `ghcr.io/quarto-dev/quarto:1.5.57` |

## Files in this workflow

```
workflows/metabolomics/
├── main.nf                       # DSL2 workflow, mode branching
├── nextflow.config                # includeConfig base.config, profiles, params.metabolomics.*
├── README.md                      # this file
└── assets/
    ├── report_template.qmd        # Quarto report template
    └── NO_FILE                    # sentinel placeholder for optional inputs

modules/metabolomics/
├── feature_table_intake.nf
├── annotation_confidence.nf
├── normalize_qc.nf
├── differential_abundance.nf
├── pathway_mapping.nf
├── standardize_output.nf
├── emit_manifest.nf
└── render_report.nf
```
