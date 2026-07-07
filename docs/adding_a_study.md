# Adding a study

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
    sample_size: 24
    biological_replicates: 3
    results_file: "data/demo/rnaseq/GIGAS_HEAT01_acute_heat_vs_control_dge_demo.tsv"
```

Each comparison needs its own `comparison_id` (unique within the study),
`tissue` and `life_stage` (from the controlled vocabularies), a
`stressor_original` (exact free text from the source — never discard this)
alongside a `stressor_standardized` ontology term, a `phenotype` term, a
`phenotype_direction`, a `resilience_classification`, and sample-size fields.
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
`registry/study_registry.csv`. Registering a `study_id` that already exists
without `--update` fails with a duplicate-study error rather than silently
overwriting.

Confirm it landed:

```bash
aree list-studies
```

## 6. Next steps

Once a study is registered, harmonize its results into the shared evidence
table — see [raw_vs_processed.md](raw_vs_processed.md) for which mode applies
and the command to run. From there it becomes eligible for
[meta-analysis](interpreting_meta_analysis.md) and
[candidate scoring](interpreting_candidate_scores.md).

## Related documentation

- [defining_a_phenotype.md](defining_a_phenotype.md)
- [resilience_vs_exposure.md](resilience_vs_exposure.md)
- [raw_vs_processed.md](raw_vs_processed.md)
- [identifier_mapping.md](identifier_mapping.md)
