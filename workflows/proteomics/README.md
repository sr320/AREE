# AREE proteomics workflow

STATUS: **both modes run end to end in CI.** `processed_results_harmonization`
runs on the demo table; `raw_reanalysis` runs for real (not stubbed) on a
simulated spectral-count matrix (`tests/fixtures/proteomics_raw`), and
`scripts/check_proteomics_raw_smoke.py` checks that the planted effects come
back and untestable proteins stay out. No container has been pulled; CI and
local runs use `-profile local`.

## What this workflow does

Implements Layer 2.C of the AREE design (`CLAUDE.md`, "Proteomics"): input
harmonization for protein abundance tables, normalization, a missingness
report, differential abundance, protein-to-gene identifier translation, and a
standardized protein-level evidence output — or, when only a processed
differential-abundance table is available, harmonizing that table directly
into the same standardized schema. Both modes end in a manifest and an HTML
report, per `docs/design.md` Sec. 7.

## The two modes

Controlled by `params.mode`, default `processed_results_harmonization`.

### `raw_reanalysis`

For studies that deposited, or published, a quantified protein x sample
matrix: intensities (MaxQuant/Skyline/DIA-NN exports) or spectral counts /
NSAF. **"Raw" here means that matrix, not mass-spec spectra (.raw/.mzML)** —
spectral search is out of scope; see `docs/roadmap.md`.

```
HARMONIZE_INPUT      wide matrix -> long; keep only the sample sheet's samples
   │
NORMALIZE            log2 on the declared input scale; zeros / values at or below
   │                 detection_limit become missing; optional median normalization
MISSINGNESS_REPORT   per-protein and per-sample missingness (never imputed)
   │
DIFFERENTIAL_ABUNDANCE  limma moderated t; a protein is tested only if detected
   │                    in >= min_detected_per_arm samples per arm; writes lfcSE
ID_TRANSLATION       gene symbol for display (UniProt or legacy CGI_ route)
   │
STANDARDIZE_OUTPUT -> EMIT_MANIFEST -> RENDER_REPORT
```

Inputs:

- `--proteomics.raw_abundance_matrix` — TSV/CSV, first column the protein
  accession, then one column per sample. Columns not in the sample sheet are
  left out of the comparison, so one matrix serves every comparison of a study.
- `--proteomics.sample_sheet` — `sample_id` and `group` (or `condition`), with
  `group` exactly `control` or `treatment` for every sample.
- `--proteomics.id_mapping_table` — defaults to the SYNTHETIC demo crosswalk;
  pass the species crosswalk (e.g. `data/reference/crosswalk/mgigas_gene_id_crosswalk.tsv`)
  for a real study. It only labels the report; evidence-level identifier
  resolution is done by `aree harmonize`.
- `--proteomics.study_id`, `--proteomics.comparison_id`.

Options (`workflows/proteomics/nextflow.config`):

| Option | Default | Meaning |
|---|---|---|
| `input_scale` | `auto` | `linear` → log2(x); `log2` → as is; `auto` keeps the old heuristic (max < 100 means log2, otherwise log2(x + 1)), which **misreads NSAF**: set `linear` for any spectral-count ratio |
| `detection_limit` | null (0 for linear) | values at or below it are missing, so a zero spectral count is "not detected" |
| `normalization` | `median` | `median` or `none`; NSAF is already normalized within each sample |
| `min_detected_per_arm` | 2 | proteins detected in fewer samples in either arm are listed in `*_detection_status.tsv` and not tested |
| `quantification` | `intensity` | `intensity` or `spectral_count`; becomes `analysis_method` (`spectral_count_limma_moderated_t`) |

### `processed_results_harmonization`

For studies that only released a processed protein-level differential
abundance table (see `registry/demo_studies/GIGAS_SAL04.yaml`):

```
STANDARDIZE_OUTPUT (validate/reshape the supplied table)
   │
EMIT_MANIFEST (warns that raw QC/normalization/missingness were not re-run)
   │
RENDER_REPORT
```

Inputs: `--proteomics.processed_results`, `--proteomics.study_id`,
`--proteomics.comparison_id`. `config/demo.config` wires this mode against the
simulated demo table. Where no missingness report exists the workflow passes the
sentinel `assets/NO_FILE`, which the modules record as "not supplied".

## Standardized output schema

`STANDARDIZE_OUTPUT` emits `<study_id>_<comparison_id>_protein_abundance_standardized.tsv`:

| Column | Notes |
|---|---|
| `protein_accession` | as in the input: UniProt accession, or a source-native ID such as an oyster_v9 `CGI_` model or ProteinProphet group |
| `gene_symbol` | display only; identity is resolved from `protein_accession` by `aree harmonize` |
| `log2FC` | treatment vs control |
| `lfcSE` | limma posterior SE in raw mode; null when a processed table does not report one (never imputed) |
| `pvalue`, `padj` | from limma, or as published |
| `missingness_percent` | percent of the comparison's samples missing the protein, before any filter |
| `analysis_method` | e.g. `spectral_count_limma_moderated_t`; null for processed tables, which then default to `MSstats` at harmonization |

## Feeding into harmonization

```bash
aree harmonize --study STUDY_ID --comparison COMPARISON_ID \
    --input path/to/STUDY_ID_COMPARISON_protein_abundance_standardized.tsv
```

`src/harmonize/proteomics.py` resolves `CGI_` IDs through the crosswalk's
`legacy_locus_tags` (`inferred`) and everything else as UniProt accessions,
carries `lfcSE` into `standard_error` and `analysis_method` per row, and adds
`processed_only` when `missingness_percent >= 10`. See
`docs/identifier_mapping.md`.

## Provenance and manifest

`EMIT_MANIFEST` writes `<study_id>_<comparison_id>_manifest.json` with the
sha256 of every input file and of the standardized table, every
`params.proteomics` value, software versions collected from each process's
`versions.yml`, the missingness summary, the limma QC (samples per arm,
proteins tested and why others were not, `df.prior`), and warnings.

## Running

```bash
# processed mode on the demo table
nextflow run workflows/proteomics -profile demo,local

# raw mode on a spectral-count (NSAF) matrix
nextflow run workflows/proteomics -profile local \
    --mode raw_reanalysis \
    --proteomics.raw_abundance_matrix path/to/nsaf_matrix.tsv \
    --proteomics.sample_sheet path/to/COMPARISON_samplesheet.tsv \
    --proteomics.input_scale linear \
    --proteomics.normalization none \
    --proteomics.quantification spectral_count \
    --proteomics.id_mapping_table data/reference/crosswalk/mgigas_gene_id_crosswalk.tsv \
    --proteomics.study_id STUDY_ID \
    --proteomics.comparison_id COMPARISON_ID
```

## What has and has not been run

- **Run in CI:** processed mode on the demo table; raw mode on the simulated
  spectral-count fixture, natively (`-profile local`), including limma and the
  Quarto report.
- **Not run:** any container image; spectral search from .raw files (out of
  scope); an intensity-scale (`input_scale log2`) real study.
