# Adding a study

> Working through several studies on different machines? Start each one with
> `aree start-study`, which claims it on its own branch. See
> [parallel_study_runs.md](parallel_study_runs.md).

This walks through registering a new public (or demo) dataset in AREE. A
"study" is registration-time metadata describing what was done, not what was
found — see [design.md](design.md#2-data-model-overview) for why studies,
comparisons, and evidence records are kept as separate tiers.

## 1. Choose single-study YAML or batch CSV

For one study, copy the template:

```bash
cp registry/studies/_TEMPLATE.yaml registry/studies/YOUR_STUDY_ID.yaml
```

For registering several studies at once, start from
`registry/studies/_batch_template.csv`, which has one row per comparison
(the same fields as the YAML template, flattened). The batch CSV is useful
for bulk entry but every study still needs to pass the same schema validation
as a YAML record — there is no separate, looser batch schema.

## 2. Fill in the required fields

The full field set is defined in `schemas/study.schema.json`. At minimum you
must provide:

- `study_id` — unique, must match the filename stem, pattern `^[A-Za-z0-9_-]+$`
- `citation`, `species`, `genome_assembly`
- `assay_type` — a list of one or more terms from
  `registry/controlled_vocabularies/assay_types.yaml` (`rnaseq`, `methylation`,
  `proteomics`, `metabolomics`, `microarray`, `qpcr`, `genotyping`)
- `analysis_mode` — `raw_reanalysis` or `processed_results_harmonization`
  (see [raw_vs_processed.md](raw_vs_processed.md))
- `comparisons` — at least one (see below)
- `data_availability` — `raw_data_available`, `processed_data_available`,
  `status` (`available` | `restricted` | `embargoed` | `unavailable`)
- `provenance` — `registered_by`, `date_registered` (ISO date)

Other useful fields: `accessions` (DOI/BioProject/GEO/SRA/ENA/ProteomeXchange/
other — all nullable), `strain_or_population`, `annotation_version`,
`platform`, `qc_status`, `analysis_status`, `quality_flags`, `limitations`.

If the study is synthetic/demo data rather than a real dataset, set
`simulated: true` and say so plainly in the `citation` and `limitations`
fields (see any file in `registry/studies/GIGAS_*.yaml` for the convention
used by the demo set — they are all prefixed `[SIMULATED]`).

## 3. The `comparisons[]` structure

A single study can report more than one phenotype or stressor — for example
`registry/studies/GIGAS_HEAT01.yaml` reports both an acute heat-shock
comparison (phenotype `survival`) and a chronic heat comparison (phenotype
`growth_under_stress`) from the same study. Rather than forcing one label per
study, phenotype/stressor/tissue/life-stage fields live at the **comparison**
level:

```yaml
comparisons:
  - comparison_id: "acute_heat_vs_control"
    tissue: "gill"
    life_stage: "adult"
    stressor_original: "Acute heat shock, 34C for 2h, versus 20C ambient seawater"
    stressor_standardized: "temperature"
    exposure_intensity: "34C"
    exposure_duration: "2h"
    exposure_timing: "acute, single exposure"
    control_condition: "20C ambient seawater, matched handling"
    treatment_condition: "34C for 2h then sampled"
    phenotype: "survival"
    phenotype_direction: "decrease"
    phenotype_unit: "percent_survival"
    resilience_classification: "resilience"
    contrast_type: "exposure_response"
    sample_size: 24
    biological_replicates: 3
    results_file: "data/demo/rnaseq/GIGAS_HEAT01_acute_heat_vs_control_dge_demo.tsv"
```

Each comparison needs its own `comparison_id` (unique within the study),
`tissue` and `life_stage` (from the controlled vocabularies), a
`stressor_original` (exact free text from the source — never discard this)
alongside a `stressor_standardized` ontology term, a `phenotype` term, a
`phenotype_direction`, a `resilience_classification`, a `contrast_type`, and
sample-size fields. `contrast_type` says whether the arms are a treatment and
its control (`exposure_response`) or groups that differ in resilience, and if
so whether they were sampled at baseline or under the challenge; see
[contrast_types.md](contrast_types.md). If a study has both kinds, choose
`meta_analysis_primary` separately within each contrast type.
`results_file` points at the processed results table for that comparison when
`analysis_mode` is `processed_results_harmonization` (or when raw-mode output
has already been produced and staged, as in the demo).

Add as many `comparisons` entries as the source study actually reports. Do
not collapse two different stressors or two different phenotypes into one
comparison — see [resilience_vs_exposure.md](resilience_vs_exposure.md) for
why this distinction matters downstream.

## 4. Validate

```bash
aree validate-study registry/studies/YOUR_STUDY_ID.yaml
```

This checks the file against `schemas/study.schema.json` and cross-checks
controlled-vocabulary fields (phenotype, stressor, tissue, life stage, assay
type, mapping/quality terms where applicable). Fix any reported errors before
proceeding; warnings are non-fatal but worth reading.

## 5. Register

```bash
aree register-study registry/studies/YOUR_STUDY_ID.yaml
```

This re-validates and appends (or, with `--update`, overwrites) a row in
`registry/study_registry.csv`. That index is generated and gitignored, so it
never appears in a PR; `aree build-registry` rebuilds it from every YAML.
Registering a `study_id` that already exists
without `--update` fails with a duplicate-study error rather than silently
overwriting.

Confirm it landed:

```bash
aree list-studies
```

## 6. Convert the published result tables

If the study is `processed_results_harmonization` mode — the common case for
public data — its results arrive as a supplementary spreadsheet or CSV that
does not match the column contract `aree harmonize` expects. Do **not** reshape
it by hand. Write an intake config instead, so the step from publication to
evidence is a command someone else can re-run and verify.

Put the downloaded artifact under `data/studies/YOUR_STUDY_ID/_source/` and
write `data/studies/YOUR_STUDY_ID/intake.yaml`:

```yaml
study_id: YOUR_STUDY_ID
assay: rnaseq

source:
  description: Supplementary Table S2 (differential expression, treated vs control)
  url: https://doi.org/...
  license: CC BY 4.0
  citation: Author A, Author B (2025). Title. Journal 1:23. https://doi.org/...
  local_copy: data/studies/YOUR_STUDY_ID/_source/supp_table_s2.xlsx
  sha256: <shasum -a 256 of that file>

output_dir: data/studies/YOUR_STUDY_ID
provenance_file: data/studies/YOUR_STUDY_ID/intake_provenance.json

conversions:
  - output_file: YOUR_STUDY_ID_treated_vs_control_dge.tsv
    source_sheet: DEGs            # spreadsheets only; omit for CSV/TSV
    column_map:                   # AREE column: source column
      gene_id: Gene
      log2FoldChange: log2FC
      lfcSE: lfcSE
      pvalue: pvalue
      padj: FDR

# Quote every note: an unquoted string containing ": " parses as a YAML mapping.
transformation_notes:
  - "Anything a reader would need to know to interpret the derived table."
```

Then:

```bash
aree intake-supplementary data/studies/YOUR_STUDY_ID/intake.yaml
```

Rules the intake step enforces, and why:

- **Map only the columns the source actually reports.** Unmapped AREE columns
  are emitted empty. If the publication reports no standard error, leave
  `lfcSE` unmapped — an empty column correctly prevents inverse-variance
  pooling, whereas a fabricated one silently corrupts every meta-analysis the
  study enters. Nothing is ever imputed.
- **Do not touch identifiers.** Decoration like `gene-LOC105331241|LOC105331241`
  is preserved verbatim; [identifier_mapping.md](identifier_mapping.md)
  describes how it is resolved later. `feature_id_original` must remain exactly
  what the authors published.
- **Do not filter beyond what the authors did.** If the published table is
  already significance-filtered, say so in `transformation_notes` — it is a
  study limitation, not something to correct here.
- **The source `sha256` is verified before every conversion.** If the publisher
  reissues the file, the intake fails rather than silently re-deriving results
  from a different version. Resolve that deliberately.

Supported sources are `.xls`, `.xlsx`, `.xlsm` (needs `pip install -e ".[intake]"`),
`.csv`, `.tsv`, and `.txt`. Verify at any later time — and in CI — with:

```bash
aree intake-supplementary data/studies/YOUR_STUDY_ID/intake.yaml --check
```

which regenerates into a temporary directory and compares checksums against the
committed files and provenance, without modifying anything.

If your source is genuinely not a differential-expression table (a methylation
region table, a protein abundance matrix, a metabolite feature table), the
intake converter does not yet cover it — see
[raw_vs_processed.md](raw_vs_processed.md) for the expected input shapes and
prepare those tables by a documented script of your own for now.

## 7. Next steps

For a `raw_reanalysis` study, `aree fetch-samplesheet` (see
[first_raw_reanalysis.md](first_raw_reanalysis.md)) writes
`fastq_manifest.tsv` alongside the sample sheet, and the two always list the
same runs:

- Each run gets exactly one manifest row per mate: 1 for `SINGLE`, 2 for
  `PAIRED`. Any other count is an error naming the runs. One exception: when
  a paired run's ENA FASTQ also includes an unpaired `RUN.fastq.gz` beside
  `RUN_1`/`RUN_2`, that file is left out and listed under
  `unpaired_files_dropped` in `ena_provenance.json`.
- URLs come from ENA's generated FASTQ (`fastq_ftp`) where it exists. Some
  public runs have none, and only the submitter's own files (`submitted_ftp`).
  When those are all FASTQ, they are used instead, with the submitter's MD5s
  and sizes. The manifest's `source` column (`fastq_ftp` or `submitted_ftp`)
  and `submitted_ftp_runs` in `ena_provenance.json` record which runs these
  are. PRJEB105019's ERR15996548 is one example. Submitted files keep their
  own names (e.g. `..._DNA_R1.fastq.gz`), so download them by
  `run_accession` and `mate`, not by file name.
- A run with neither (no ENA FASTQ, and submitted files that are BAM, CRAM
  or missing) is an error naming the run, not a silent omission. Resolve it
  with the submitter or ENA before going further. Do not hand-edit the
  manifest.

**Before downloading** a `raw_reanalysis` study, pre-screen its replicates on
the first ~1 M reads of each library:

```bash
python scripts/prescreen_replicates.py --source ncbi \
    --samplesheet data/studies/STUDY_ID/samplesheet.tsv \
    --manifest data/studies/STUDY_ID/fastq_manifest.tsv \
    --salmon-index data/reference/GCF_963853765.1/salmon_index \
    --level CONTROL_LEVEL --level TREATMENT_LEVEL
```

It quantifies only those reads and applies the same dispersion statistic as
the full check below, at about 1 GB and 10 minutes per study instead of tens
of GB. `--source ncbi` uses `fastq-dump -X` (sra-tools, in
`environment.yml`) and is usually much faster than ENA. Use `--source ena`
to screen the exact files the manifest would download. It was calibrated on
two studies with known outcomes (1 M reads per library; full-depth values in
brackets):

| Study | Pre-screen | Full depth | Verdict |
|---|---|---|---|
| `WANG2024_HEAT_GILL` | 0.0038 / 0.0029 | (0.0030 / 0.0009) | fail |
| `LI2022_HEAT_PLOIDY` d_0 / d_f | 0.0191 / 0.0218 | (0.0139 / 0.0143) | pass |

The subsample reads a little high, so a borderline study can pass the
pre-screen and still fail at full depth. A pass is permission to download,
not a substitute for the full check. A fail means: do not download the study
as is, and record why in [candidate_studies.md](candidate_studies.md).

Check after the workflow runs and **before
harmonizing** that its replicates are independent animals:

```bash
python scripts/check_replicate_dispersion.py \
    --samplesheet data/studies/STUDY_ID/samplesheet.tsv \
    --salmon-dir RESULTS_DIR/COMPARISON_ID/rnaseq/salmon
```

It exits non-zero if any group's median gene dispersion falls below 0.005.
Real *M. gigas* biological replicates sit around 0.03–0.05. A group far below
that is behaving like technical replicates, and its DESeq2 standard errors
are too small to pool. `IOCAS2022_OA_ENERGY` and `WANG2024_HEAT_GILL` failed
this check even though the deposited metadata labelled every library a
biological replicate; see [candidate_studies.md](candidate_studies.md). If a
study fails, set `qc_status: failed`, record the evidence, and do not
harmonize it. For individually sequenced animals, comparing allele fractions
at expressed SNPs can show *why* a group failed (one RNA source vs distinct
animals). `data/studies/WANG2024_HEAT_GILL/genotype_check/` is a worked
example. It cannot tell pools of larvae from one cross apart, since those
share parental genotypes.

Once a study is registered and has passed QC, harmonize its results into the shared evidence
table — see [raw_vs_processed.md](raw_vs_processed.md) for which mode applies
and the command to run. From there it becomes eligible for
[meta-analysis](interpreting_meta_analysis.md) and
[candidate scoring](interpreting_candidate_scores.md).

## Related documentation

- [defining_a_phenotype.md](defining_a_phenotype.md)
- [resilience_vs_exposure.md](resilience_vs_exposure.md)
- [raw_vs_processed.md](raw_vs_processed.md)
- [identifier_mapping.md](identifier_mapping.md)
