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
  Python. Only records with a usable standard error contribute to the pooled
  estimate or replication count; lower-confidence aliases within one
  comparison are de-duplicated, and multiple contrasts from one study fail
  closed until a covariance-aware model or prespecified contrast is supplied.
  `aree meta-analyze` runs against the demo evidence table and produces pooled
  results, including the conflicting-evidence case documented in
  [interpreting_meta_analysis.md](interpreting_meta_analysis.md).
- **Candidate prioritization.** `src/prioritize/scoring.py` and
  `src/prioritize/rank.py` implement the full transparent scoring formula and
  the three-tier hard-gated ranking system. Simulation status and species taxid
  remain partition keys through ranking and reporting; `aree build-evidence-cards`
  generates a markdown card per candidate plus an `index.json` under
  `reports/evidence_cards/`.
- **CLI.** All six commands documented in the README/design docs
  (`validate-study`, `register-study`, `list-studies`, `harmonize`,
  `meta-analyze`, `build-evidence-cards`) are implemented in
  `src/aree/cli.py` and run against the demo data on a clean install.
- **Species and assembly reference data.**
  `registry/controlled_vocabularies/species.yaml` maps scientific names and
  accepted synonyms (*Crassostrea*/*Magallana gigas*) onto NCBI taxids, and
  `data/reference/genome_assemblies.yaml` carries NCBI-verified accessions for
  both oyster assemblies in use. Both are enforced by `aree validate-study`,
  and evidence records now carry `species_taxid` plus
  `identifier_annotation_release` so that the crossing between a study's
  assembly and the annotation its identifiers were resolved against is
  visible — see
  [handling_genome_versions.md](handling_genome_versions.md).
- **Automated tests.** `tests/` covers schema validation, registry
  (duplicate-ID handling), identifier mapping, meta-analysis, prioritization,
  evidence-card generation, origin/species collision regression cases, and a
  full demo-pipeline integration test.
- **CI, interface, and documentation.** `.github/workflows/ci.yml` runs Python,
  real-study, schema, documentation, and Nextflow smoke jobs; `app/main.py`
  provides the Streamlit browser; `docs/_quarto.yml` renders this site.

## Scaffolded but not yet production-ready

- **Nextflow workflows.** All four now run `processed_results_harmonization`
  end to end, and `workflows/rnaseq` has additionally run the full
  `raw_reanalysis` path against real public FASTQ (PRJNA1329250, subsampled;
  see [first_raw_reanalysis.md](first_raw_reanalysis.md)). None of this was
  true before 2026-08-28 — three of the four did not compile on current
  Nextflow. Still outstanding: full-depth runs, the raw paths for methylation,
  proteomics and metabolomics, and any execution inside the declared
  containers. CI executes all four processed paths on Nextflow 26.04 and uses
  tiny paired FASTQ fixtures plus `-stub-run` to exercise the RNA-seq raw DAG;
  the stub run verifies wiring and artifact contracts, not tool execution or
  biological validity.
- **Containers.** `containers/README.md` documents the intended
  Docker/Apptainer approach, but container images/definitions themselves are
  not yet built.
- **Ortholog mapping.** A real 33,356-gene *M. gigas* identifier crosswalk is
  committed and used for real studies, but its `orthogroup_id` field remains
  empty. Cross-species pooling still needs a real OrthoFinder/OrthoDB-derived
  orthology layer — see [identifier_mapping.md](identifier_mapping.md).

## Next milestones

- **First full-depth real raw reanalysis.** Run all prespecified
  `CALLA2026_OSHV` contrasts at full depth, but do not combine its repeated
  contrasts as independent effects. Select one comparison per study/phenotype
  or implement a covariance-aware within-study model first.
- **First real cross-study pool.** Curate and reanalyze a compatible independent
  pathogen-challenge study such as `PRJNA593309`, aligning phenotype and
  contrast definitions before pooling.
- **Container verification.** Execute the CI fixtures with the declared Docker
  images, pin an AREE utility image, and record introspected rather than declared
  software versions.
- **Additional species.** The schema and vocabularies are species-aware and
  ranking/reporting partitions on `species_taxid`; onboarding another species
  still requires its reference data and identifier crosswalk.
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
