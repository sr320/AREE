# Technical architecture

AREE is organized as five connected layers, following [design.md](design.md).
This page describes how they connect and which directories implement each
one.

```
Layer 1: Study Registry & Intake
        |  (registry/, src/intake)
        v
Layer 2: Standardized Reanalysis Workflows
        |  (workflows/, modules/, containers/)
        v
Layer 3: Cross-study Harmonization
        |  (src/harmonize)
        v
Layer 4: Meta-analysis & Prioritization
        |  (src/meta_analysis, src/prioritize)
        v
Layer 5: User-facing Outputs
           (app/, docs/, reports/)
```

## Layer 1: Study Registry and Dataset Intake

**Directories:** `registry/studies/` (one YAML per study, plus
`_TEMPLATE.yaml` and `_batch_template.csv`), `registry/controlled_vocabularies/`
(phenotype, stressor, assay-type, tissue, life-stage, mapping-confidence, and
quality-flag ontologies), `registry/study_registry.csv` (flat index of
registered studies), `schemas/study.schema.json` (validation contract),
`src/intake/` (`schema_validate.py`, `registry.py`).

A study is registered by validating a YAML file against the JSON Schema and
controlled vocabularies (`aree validate-study`), then appending it to
`registry/study_registry.csv` (`aree register-study`). This layer records
facts about what was done — species, genome assembly, assay type, treatment
vs. control, sample sizes, data availability — never derived results. See
[adding_a_study.md](adding_a_study.md).

## Layer 2: Standardized Reanalysis Workflows

**Directories:** `workflows/{rnaseq,methylation,proteomics,metabolomics}/`,
`modules/`, `containers/`, `config/`.

Nextflow workflows for raw-data reanalysis: FASTQ/read QC, alignment or
pseudoalignment, quantification or methylation calling, differential
analysis, and standardized result-table emission, per assay type. Every
workflow emits a machine-readable manifest (parameters, tool versions, input
checksums, QC metrics) alongside its result table. The RNA-seq raw path has
executed end to end against subsampled real FASTQ, and all four processed-result
paths execute end to end. Current-version CI runs those processed paths and
structurally exercises the RNA-seq raw DAG with tiny fixtures and explicit
stubs. Container profiles and the other raw paths remain unverified; see
[roadmap.md](roadmap.md).

## Layer 3: Cross-study Harmonization

**Directories:** `src/harmonize/` (`core.py`, `identifiers.py`, `schema.py`),
`data/mappings/` (`gene_id_crosswalk.tsv`, `ambiguous_symbol_map.yaml`),
`schemas/evidence.schema.json`.

This is where assay-specific results (a DE table, a DMR table, a protein
abundance table, a metabolite feature table) are converted into the shared,
assay-agnostic evidence schema — one row per molecular feature per
comparison. `src/harmonize/identifiers.py` resolves feature identifiers
through the documented hierarchy (NCBI Gene ID > Ensembl > UniProt > locus ID
> gene symbol > orthogroup) and assigns a `mapping_confidence`. Output is
appended to `reports/evidence/evidence_table.tsv`, invoked via
`aree harmonize --study STUDY_ID [--input path/to/results.tsv]`. See
[identifier_mapping.md](identifier_mapping.md) and
[raw_vs_processed.md](raw_vs_processed.md).

## Layer 4: Meta-analysis and Candidate Prioritization

**Directories:** `src/meta_analysis/` (`pooling.py`, `effect_sizes.py`,
`run.py`), `src/prioritize/` (`scoring.py`, `rank.py`), `src/reporting/`
(`evidence_cards.py`).

`aree meta-analyze` groups the evidence table by standardized feature,
phenotype, and feature type, and pools effect sizes with a DerSimonian-Laird
random-effects estimator (`src/meta_analysis/pooling.py`), computed directly
in Python rather than via an R meta-analysis package (see
[methods.md](methods.md)). Output feeds `src/prioritize`, which computes a
transparent, weighted candidate score (`scoring.py`) and assigns one of three
tiers using hard gates that a score cannot override (`rank.py`). See
[interpreting_meta_analysis.md](interpreting_meta_analysis.md) and
[interpreting_candidate_scores.md](interpreting_candidate_scores.md).
`aree build-evidence-cards` (`src/reporting/evidence_cards.py`) ranks every
candidate into `reports/evidence_cards/candidates.tsv` and renders a markdown
evidence card for each candidate with a significant signal (pooled or
study-level BH-adjusted p ≤ 0.05 by default; `--all-cards` renders all), plus
an `index.json` of the cards written. Ranking and card rendering are a single
indexed pass over the evidence table, so a genome-wide pool of ~30,000
candidates ranks in seconds.

## Layer 5: User-facing Outputs

**Directories:** `app/` (Streamlit interface), `docs/` (this documentation
set, renderable as a Quarto site), `reports/` (generated evidence table,
meta-analysis tables, evidence cards — build artifacts, not source data).

The interface layer reads only already-generated files under `registry/` and
`reports/` — it does not perform analysis itself. `streamlit run app/main.py`
launches the interactive filter/search layer; `quarto render docs/` builds
the static documentation/report site. Neither requires login or a hosted
deployment; both operate on the local filesystem. See
[installation.md](installation.md) for how to launch each.

## Cross-cutting concerns

- **Provenance** (source accession, checksums, tool versions, workflow
  version, `date_generated`, `generated_by`) is threaded through every layer
  from intake to evidence card — see
  [governance_and_provenance.md](governance_and_provenance.md).
- **Controlled vocabularies** (`registry/controlled_vocabularies/`) are
  shared across Layer 1 (study registration) and Layer 3 (evidence records),
  so a phenotype or stressor term means the same thing at both stages.
- **`common`** (`src/common/`) centralizes shared paths (`REPO_ROOT`,
  `REGISTRY_DIR`, `EVIDENCE_TABLE_PATH`, etc.), YAML/JSON IO, checksum
  helpers, and controlled-vocabulary loaders used by every other package.

## Related documentation

- [design.md](design.md) — full data model and rationale behind this layering
- [roadmap.md](roadmap.md) — what is complete/runnable vs. scaffolded per layer
