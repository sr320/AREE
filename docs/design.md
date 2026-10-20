# AREE Design Document

Status: Phase 1 design, MVP-scoped. This document fixes the data model, controlled
vocabularies, provenance model, candidate-ranking framework, and repository
architecture that the rest of the codebase implements. Where the MVP takes a
shortcut relative to the long-term vision in the proposal, that is called out
explicitly under **Assumptions**.

## 1. Problem framing

Public *Crassostrea gigas* (and future shellfish/aquaculture) omics datasets are
scattered across repositories, genome versions, assay types, and inconsistent
phenotype/stressor vocabularies. A statistically significant hit in one study is
not evidence of a validated biomarker. AREE's job is to:

1. make heterogeneous studies comparable without pretending they are identical,
2. keep every transformation traceable back to its source,
3. surface convergence (or contradiction) across studies and molecular layers,
4. and refuse to silently upgrade "suggestive" into "validated."

AREE is an **evidence engine**, not a claims engine. Every output layer down to
the evidence card must let a skeptical reviewer trace a number back to a source
file and a set of parameters.

## 2. Data model overview {#2-data-model-overview}

Three tiers, each with its own identity and lifecycle:

```
Study (registry/studies/*.yaml)
  └── one or more Comparisons (treatment vs control, within a study)
        └── many Evidence Records (one row per feature per comparison)
              └── rolled up into Candidates (cross-study synthesis)
```

- **Study**: registration-time metadata. Describes what was done, not what was
  found. Immutable once registered except for curation corrections (tracked).
- **Evidence record**: the atomic unit of the harmonized evidence table — one
  molecular feature's effect estimate in one comparison from one study.
  Assay-agnostic schema; assay-specific detail is preserved in `source_file` and
  workflow manifests, not flattened away.
- **Candidate**: a synthesized view across evidence records sharing a
  standardized feature identifier (or ortholog group), computed by the
  meta-analysis / prioritization layer. Candidates are derived, not curated —
  they must be reproducible from the evidence table alone.

### Why this shape

Phenotypes, stressors, tissues, and life stages are properties of the
**comparison**, not the study (a single study can report multiple phenotypes,
e.g. survival AND growth, under multiple stressors). The schema therefore
carries these fields at the evidence-record level (per comparison), while the
study record carries the union/description of what the study covers. This
avoids collapsing multi-phenotype or multi-stressor studies into one label.

## 3. Controlled vocabularies {#3-controlled-vocabularies}

Defined in `registry/controlled_vocabularies/`:

- `phenotype_ontology.yaml` — resilience-related phenotype terms (survival,
  thermal tolerance, disease resistance, etc.), each tagged with a
  `resilience_relevance` class (`resilience`, `stress_response`, `disease`,
  `exposure_only`) so the system never conflates "gene moved under stress" with
  "gene predicts resilience."
- `stressor_ontology.yaml` — environmental/experimental stressor classes.
- `assay_types.yaml`, `feature_types.yaml`, `tissue_types.yaml`,
  `life_stages.yaml`, `mapping_confidence.yaml`, `quality_flags.yaml` — smaller
  enumerations reused across schemas.

Original free-text values (e.g. the exact treatment description from a paper)
are always preserved alongside the controlled-vocabulary mapping — see
`treatment_original` vs `stressor_standardized` in the study schema. Mapping is
additive, never destructive.

## 4. Provenance model {#4-provenance-model}

Every workflow output and every evidence record carries:

- source accession(s) and study_id
- input file name(s) and checksum(s) (sha256)
- parameter set (explicit key/value, not "defaults")
- workflow name + semantic version
- tool versions (from container manifest or environment lock)
- reference genome/annotation version
- `date_generated` (ISO 8601, supplied by caller — the codebase does not call
  wall-clock time internally, so provenance timestamps are always explicit and
  reproducible)
- `generated_by` (`automated:<workflow_id>@<version>` or `curator:<name>`)
- manual curation decisions, if any, as a separate append-only log entry

Missingness is a first-class value: fields the source study does not report are
recorded as `null` with a `missing_reason` where feasible, never silently
dropped from the schema or imputed.

## 5. Identifier harmonization {#5-identifier-harmonization}

Identifier hierarchy (highest to lowest precedence when multiple are available
for the same feature): **NCBI Gene ID → Ensembl gene ID → UniProt accession →
locus ID from the reference oyster annotation (e.g. LOC/CGI) → gene symbol →
orthogroup**. The *original* identifier as reported by the source study is
always retained in `feature_id_original`; `feature_id_standardized` is the
harmonized value; `mapping_confidence` records how much to trust the
translation:

`exact | one_to_one_ortholog | one_to_many_ortholog | many_to_one_ortholog | inferred | unresolved`

`unresolved` is a valid, expected outcome — it is not an error state, and
downstream meta-analysis code must handle it (by excluding the record from
identifier-level pooling while still keeping it in the evidence table).

## 6. Candidate scoring framework (transparent, not black-box) {#6-candidate-scoring-framework-transparent-not-black-box}

Score is a weighted sum of named, independently inspectable components, all
computed from fields already in the evidence/candidate tables:

| Component | What it measures | Source field(s) |
|---|---|---|
| `n_studies_score` | independent studies supporting the candidate | count of distinct `study_id` |
| `sample_size_score` | total biological replication behind the pooled estimate | sum of `sample_size` |
| `effect_magnitude_score` | typical absolute standardized effect size | pooled `effect_size` |
| `significance_score` | adjusted-significance strength | pooled/median `adjusted_p_value` |
| `direction_consistency_score` | agreement in direction of effect across studies | fraction of records matching majority direction |
| `phenotype_relevance_score` | resilience vs. stress-response vs. exposure-only weighting | most conservative contributing `resilience_classification` |
| `context_breadth_score` | spread across tissues/life stages | distinct `tissue` × `life_stage` combinations |
| `assay_diversity_score` | molecular layers with significant same-phenotype support for the feature | `feature_type` → `molecular_layer` (feature_types vocabulary), `adjusted_p_value` |
| `mapping_confidence_score` | how trustworthy the identifier harmonization is | `mapping_confidence` |
| `quality_score` | study/data quality flags | `quality_flags` |
| `heterogeneity_penalty` | subtracted; large cross-study inconsistency (I²) reduces score | meta-analysis `I2` |

Weights and the exact formula live in `src/prioritize/scoring.py` as named
constants with inline rationale — this is intentionally the single place to
audit or tune the score. The function is pure (same inputs → same score) so
`candidate-score reproducibility` is directly testable.

### Ranking tiers

1. **High-priority cross-study candidates** — ≥2 independent studies,
   interpretable phenotype (not exposure-only), direction-consistency ≥ 0.7,
   acceptable quality flags, and a pooled effect that survives
   Benjamini–Hochberg control within its phenotype/feature-type family
   (adjusted p ≤ 0.05).
2. **Multi-omics convergence candidates** — ≥2 distinct molecular layers
   (as declared per feature type in the `feature_types` vocabulary) each with
   a significant record for the same standardized feature/orthogroup under the
   same phenotype and contrast type, the candidate's own layer among them, with the mapping path
   shown explicitly in the evidence card (this is why
   `feature_id_standardized` and `mapping_confidence` are mandatory, not
   optional, fields). Cross-phenotype overlap is shown on the card but not
   credited.
3. **Metabolite-linked convergence candidates** — significant same-phenotype,
   same-contrast evidence connected by KEGG reaction annotation or a curated
   transporter–substrate link. A separate tier below shared-identity convergence,
   with no contribution to the numeric score; currency and hub compounds are
   excluded. See [Metabolite–gene links](metabolite_gene_links.md).
4. **Emerging candidates** — single-study support, but biologically plausible
   or a strong effect; always explicitly labeled "requires replication" and
   assigned when none of the stronger evidence gates is met, regardless of score.

A candidate never moves from "emerging" to a higher tier automatically just by
score — the study-count and assay-diversity gates above are hard requirements,
not just scoring inputs. This is the mechanism that prevents "significant in
one study" from reading as "validated."

## 7. Raw vs. processed-results modes {#7-raw-vs-processed-results-modes}

Every study registration declares `analysis_mode`: `raw_reanalysis` or
`processed_results_harmonization`. Raw mode expects the Nextflow workflows in
`workflows/` to run against FASTQ/BAM/etc. and emit the standardized result
schema. Processed mode expects a results file (DE table, DMR table, protein
abundance table, feature table) that `src/harmonize` maps directly into
the same evidence schema, with a `quality_flags` entry noting that raw QC could
not be verified independently. Both modes converge on the identical evidence
schema — meta-analysis code does not need to know which mode produced a row.

## 8. Repository architecture

See repository tree in [README.md](https://github.com/sr320/AREE/blob/main/README.md). Key separation:

- `registry/` — data *about* studies (facts, never derived).
- `workflows/` + `modules/` — Nextflow workflows for raw-data reanalysis.
- `src/` — the Python packages: CLI (`aree`), intake/validation, harmonization,
  meta-analysis, prioritization, evidence-card reporting. This is where
  "comparable evidence" is actually produced.
- `data/demo/` — synthetic, clearly labeled demo inputs so the whole pipeline
  runs without any real download.
- `app/` — Streamlit interface over the harmonized outputs.
- `docs/` — Quarto site + narrative documentation.

## 9. Explicit assumptions (MVP scope) {#9-explicit-assumptions-mvp-scope}

- Demo data is entirely synthetic and is labeled `"simulated": true` in every
  demo study record and in filenames (`*_demo.tsv`). No real accession numbers
  are invented; the final response lists *real* candidate datasets to curate
  next, described by type/context rather than fabricated IDs.
- RNA-seq has full-depth real runs; methylation has a real Bismark/pooled-limma
  run, and proteomics has a real quantified-matrix run. CI also exercises
  metabolomics fixtures. Container execution remains unverified; see
  [implementation status](implementation_status.md).
- Demo identifier crosswalks are synthetic. Real studies use species-specific
  crosswalks and a versioned OrthoFinder run; cross-species pooling is separate
  from within-species gene pooling.
- Meta-analysis uses a standard random-effects (DerSimonian–Laird) estimator
  implemented directly in Python for auditability, rather than depending on a
  compiled R meta package, so `aree meta-analyze` has no R runtime dependency
  in the MVP. R is still used for the optional DESeq2/edgeR-style workflow
  scaffolds under `workflows/rnaseq/`.
- The demo covers oyster and mussel. Real evidence currently covers oyster
  and Manila clam. Explicit species and assembly fields preserve partitions
  (see [Adding a species](adding_a_species.md)).
- The Streamlit app and Quarto docs both ship; Quarto is the primary
  documentation/report site, Streamlit is the interactive filter/search layer.
- No authenticated hosted application is provided. The public Quarto site
  is deployed to GitHub Pages; the CLI and Streamlit app read local files.
