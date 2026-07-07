# Roadmap

This page distinguishes what is built and runnable today from what remains
scaffolded or planned. It is written to be checked against the repository,
not aspirational — see [design.md](design.md#9-explicit-assumptions-mvp-scope)
for the assumptions this roadmap is consistent with.

## MVP: complete and runnable today

- **Schemas.** `schemas/study.schema.json` and `schemas/evidence.schema.json`
  are complete JSON Schema definitions, enforced by `aree validate-study` and
  by `src/harmonize`.
- **Controlled vocabularies.** Phenotype ontology (15 terms), stressor
  ontology (11 terms), assay types, tissue types, life stages, mapping
  confidence, and quality flags all exist under
  `registry/controlled_vocabularies/` and are actively checked/consumed by
  code (`src/common/load_vocab`, `src/intake/schema_validate.py`).
- **Study registry.** Six demo studies are registered
  (`registry/studies/GIGAS_{HEAT01,OA02,PATH03,SAL04,LARV05,GROW06}.yaml`),
  covering the required demo mix: two RNA-seq raw-reanalysis studies
  (`GIGAS_HEAT01`, `GIGAS_OA02`), one methylation study (`GIGAS_PATH03`), one
  proteomics processed-results-only study (`GIGAS_SAL04`), one metabolomics
  processed-results-only study (`GIGAS_GROW06`), one imperfect-identifier
  example (`GIGAS_LARV05`, legacy symbol-only annotation), and one
  conflicting-direction example (`sod1`/`LOC105331241` across `GIGAS_OA02`
  and `GIGAS_LARV05`). `registry/study_registry.csv` is the resulting index.
  `aree validate-study` and `aree register-study` both work end to end
  against these files.
- **Harmonization for all four assay types.** `src/harmonize` converts
  processed RNA-seq, methylation, proteomics, and metabolomics result tables
  into the shared evidence schema; `aree harmonize` runs successfully against
  every demo study and writes to `reports/evidence/evidence_table.tsv`.
- **Identifier harmonization.** `src/harmonize/identifiers.py` implements the
  full precedence hierarchy and all six `mapping_confidence` levels, backed
  by a synthetic but internally consistent crosswalk
  (`data/mappings/gene_id_crosswalk.tsv`) and an ambiguous-symbol exception
  table (`data/mappings/ambiguous_symbol_map.yaml`).
- **Meta-analysis.** `src/meta_analysis/pooling.py` implements
  DerSimonian-Laird random-effects pooling with Q, I², and tau² directly in
  Python; `aree meta-analyze` runs against the demo evidence table and
  produces real pooled results, including the genuine conflicting-evidence
  case documented in
  [interpreting_meta_analysis.md](interpreting_meta_analysis.md).
- **Candidate prioritization.** `src/prioritize/scoring.py` and
  `src/prioritize/rank.py` implement the full transparent scoring formula and
  the three-tier hard-gated ranking system; `aree build-evidence-cards`
  generates a markdown card per candidate plus an `index.json` under
  `reports/evidence_cards/`.
- **CLI.** All six commands documented in the README/design docs
  (`validate-study`, `register-study`, `list-studies`, `harmonize`,
  `meta-analyze`, `build-evidence-cards`) are implemented in
  `src/aree/cli.py` and run against the demo data on a clean install.
- **Automated tests.** `tests/` covers schema validation, registry
  (duplicate-ID handling), identifier mapping, meta-analysis, prioritization,
  evidence-card generation, and a full demo-pipeline integration test.

## Scaffolded but not yet production-ready

- **Nextflow raw-data workflows.** `workflows/{rnaseq,methylation,proteomics,metabolomics}/`
  each contain a `main.nf`, `nextflow.config`, and a Quarto report template
  (`assets/report_template.qmd`). These are structurally complete scaffolds
  — channel definitions, process stubs, config, container references — but
  have **not been executed against real FASTQ/BAM/spectra** in this build.
  The demo's `raw_reanalysis`-mode studies (`GIGAS_HEAT01`, `GIGAS_OA02`,
  `GIGAS_PATH03`) use pre-staged standardized result tables under
  `data/demo/` rather than live Nextflow execution — there was no raw data or
  compute budget available to run them for real in this environment. Treat
  these as a documented starting point for real execution, not a validated
  pipeline.
- **Containers.** `containers/README.md` documents the intended
  Docker/Apptainer approach, but container images/definitions themselves are
  not yet built.
- **Ortholog mapping.** The identifier crosswalk is a small, illustrative,
  entirely synthetic table (21 genes), not a full OrthoFinder/OrthoDB-based
  ortholog-calling pipeline. It is sufficient to demonstrate the
  mapping-confidence mechanism but is not a real reference crosswalk — see
  [identifier_mapping.md](identifier_mapping.md).
- **Genome assembly reference.** `data/reference/genome_assemblies.yaml`
  currently holds one placeholder entry with an explicitly fake
  `ncbi_assembly_accession` (`PLACEHOLDER_CONFIRM_BEFORE_REAL_USE`) — see
  [handling_genome_versions.md](handling_genome_versions.md).

## Planned / not yet started

- **CI.** No `.github/workflows/` are present in this build yet. The design
  calls for GitHub Actions covering Python tests, R checks, schema
  validation, demo-report rendering, linting, and end-to-end demo-command
  checks — these are specified but not implemented.
- **Streamlit app.** `app/` does not yet contain a `main.py`; the interactive
  filter/search interface described in [architecture.md](architecture.md#layer-5-user-facing-outputs)
  is planned but not built in this pass. The CLI and the generated
  `reports/` artifacts (evidence table, meta-analysis tables, evidence cards)
  are usable in the meantime without the interactive app.
- **Quarto site configuration.** This `docs/` directory contains the
  narrative and reference Markdown content, but `_quarto.yml`/site
  configuration for rendering it as a browsable static site is handled
  separately and is not part of this documentation pass.
- **Real public datasets.** All six registered studies are synthetic demo
  data, clearly labeled `simulated: true`. No real public *C. gigas* dataset
  has been curated into the registry yet.
- **Additional species.** The schema and vocabularies are species-agnostic
  (see [adding_a_species.md](adding_a_species.md)), but onboarding a second
  species also requires a small code change to support multiple identifier
  crosswalks (`src/harmonize/identifiers.py` currently points at one fixed
  crosswalk path) — this has not been implemented.
- **Expanded ontologies.** Phenotype, stressor, tissue, and life-stage
  vocabularies cover the terms specified in the founding proposal but will
  need extension as more species and study types are added (e.g. finfish
  tissue types, additional stressor classes).
- **Hosted/authenticated deployment.** Explicitly out of scope for the first
  release per [design.md](design.md#9-explicit-assumptions-mvp-scope) — no
  login or cloud dependency is planned for the near term.

## Related documentation

- [design.md](design.md) — assumptions this roadmap tracks against
- [architecture.md](architecture.md) — which directories implement which layer
