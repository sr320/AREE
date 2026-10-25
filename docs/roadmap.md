# Roadmap

This page distinguishes what is built and runnable today from what remains
scaffolded or planned. It is written to be checked against the repository,
not aspirational — see [design.md](design.md#9-explicit-assumptions-mvp-scope)
for the assumptions this roadmap is consistent with.

## MVP: complete and runnable today

- **Schemas.** `schemas/study.schema.json` and `schemas/evidence.schema.json`
  are complete JSON Schema definitions, enforced by `aree validate-study` and
  by `src/harmonize`.
- **Controlled vocabularies.** Phenotype ontology (18 terms), stressor
  ontology (12 terms), assay types, tissue types, life stages, mapping
  confidence, and quality flags all exist under
  `registry/controlled_vocabularies/` and are actively checked/consumed by
  code (`src/common/load_vocab`, `src/intake/schema_validate.py`).
- **Study registry.** Seven demo studies are registered: six oyster
  (`registry/demo_studies/GIGAS_{HEAT01,OA02,PATH03,SAL04,LARV05,GROW06}.yaml`) and
  one mussel (`MYTGAL_HEAT01`) that exercises species partitioning. The six
  oyster studies cover the required demo mix: two RNA-seq raw-reanalysis studies
  (`GIGAS_HEAT01`, `GIGAS_OA02`), one methylation study (`GIGAS_PATH03`), one
  proteomics processed-results-only study (`GIGAS_SAL04`), one metabolomics
  processed-results-only study (`GIGAS_GROW06`), one imperfect-identifier
  example (`GIGAS_LARV05`, legacy symbol-only annotation), and one
  conflicting-direction example (`sod1`/`LOC105331241` across `GIGAS_OA02`
  and `GIGAS_LARV05`). They live in `registry/demo_studies/`, apart from the real studies, and
  index into `registry/demo_study_registry.csv`; `registry/study_registry.csv`
  lists real studies only. Both are rebuilt by `aree build-registry` rather
  than committed.
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
  the four-tier hard-gated ranking system. Simulation status and species taxid
  remain partition keys through ranking and reporting; the top tier additionally
  requires a BH-adjusted pooled p ≤ 0.05. `aree build-evidence-cards` ranks every
  candidate into `reports/evidence_cards/candidates.tsv` and renders a markdown
  card (plus `index.json`) for each candidate with a significant signal.
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
- **Query layer, Explore page and MCP server.** `src/query/` is one read-only
  query layer over the evidence, candidates, cards and registry, read from
  `reports/` or from the Parquet tables the site publishes. The Explore page
  queries those tables in the browser, and `aree mcp` exposes the same queries
  as nine tools to Claude and other MCP clients, locally over stdio, from a
  pipeline run or from the published tables (`aree fetch-site-data`). See
  [querying_aree.md](querying_aree.md).

## Scaffolded but not yet production-ready

- **Nextflow workflows.** All four now run `processed_results_harmonization`
  end to end, and `workflows/rnaseq` has additionally run the full
  `raw_reanalysis` path at full depth against real public FASTQ (including PRJNA1329250;
  see [first_raw_reanalysis.md](first_raw_reanalysis.md)). None of this was
  true before 2026-08-28 — three of the four did not compile on current
  Nextflow. Methylation now has a real CIBNOR2021 run using Bismark and
  pooled family-level limma analysis; proteomics has a real quantified-matrix
  run for TIMMINS2014. Metabolomics has fixture execution but no real study
  that passed QC. Execution inside the declared containers remains unverified. CI executes all four processed paths on Nextflow 26.04 and uses
  tiny paired FASTQ fixtures plus `-stub-run` to exercise the RNA-seq raw DAG;
  the stub run verifies wiring and artifact contracts, not tool execution or
  biological validity.
- **Containers.** `containers/README.md` documents the intended
  Docker/Apptainer approach, but container images/definitions themselves are
  not yet built.
- **Ortholog mapping and cross-species evidence.** Orthogroups across the seven
  RefSeq-annotated species come from an OrthoFinder run and are attached to
  evidence at harmonization. Cross-species meta-analysis, cross-species
  candidates with their cards, and the dashboard gate exist; on real data the
  oyster and Manila clam pool through 5,975 single-copy orthogroups on thermal
  tolerance, none convergent because the clam studies are exposure-only. See
  [identifier_mapping.md](identifier_mapping.md#orthogroups).
- **Species beyond the oyster.** Seven mussel and clam species are registered
  with crosswalks: six from NCBI Gene and one, *M. coruscus*, from its
  authors' GFF with a recorded homology search. The Manila clam has the first
  real within-species pool outside the oyster (k = 2). The dashboard and app
  show per-species counts, species filters and species identity
  ([species_identity.md](species_identity.md)).
- **Ask AREE public chat.** `aree chat-serve` (`src/chat/`) lets Claude answer
  website questions by calling the query layer's nine tools, on real studies
  only, behind a daily spend budget, per-visitor limits and an optional
  Cloudflare Turnstile bot check. It is packaged as a container, built and
  health-checked in CI, and redeployed to Cloud Run after each site build. On
  2026-10-25 the deployed service is healthy but paused (its daily budget
  variable is 0), and the bot check is not yet on. Its answer-quality eval
  passes 69% of answers on every check, the failures all on faithfulness. See
  [ask_aree.md](ask_aree.md).

## Next milestones

- **Remaining raw reanalyses.** Complete the unrun prespecified
  `CALLA2026_OSHV` contrasts at full depth, but do not combine its repeated
  contrasts as independent effects. Select one comparison per study/phenotype
  or implement a covariance-aware within-study model first.
- **Broader independent replication.** Real cross-study pools already exist
  for oyster disease susceptibility and oyster and Manila clam thermal
  tolerance. Add compatible studies with measured outcomes, preserving
  species and contrast-type partitions; see the [dashboard](index.qmd).
- **Container verification.** Execute the CI fixtures with the declared Docker
  images and pin an AREE utility image. Current RNA-seq runs collect observed
  process versions; older manifests still need the provenance caveat documented
  on the implementation-status page.
- **First real mussel studies.** No real mussel study is harmonized yet; the
  WP5 backlog's *M. edulis* heat pair is the nearest within-species pool.
- **A clam study that measured tolerance**, so a cross-species row can reach
  convergence on a tolerance phenotype rather than on exposure alone.
- **GFF-route species in an orthogroup run**, and crosswalks for razor clam
  and geoduck once their authors' annotation files are sourced.
- **Expanded ontologies.** Phenotype, stressor, tissue, and life-stage
  vocabularies cover the terms specified in the founding proposal but will
  need extension as more species and study types are added (e.g. finfish
  tissue types, additional stressor classes).
- **Reopen Ask AREE safely.** Turn on the Turnstile bot check, confirm in
  Cloud Run logs that the last `X-Forwarded-For` entry is the visitor's
  address, then raise the daily budget above 0. Before further prompt tuning,
  grow the eval beyond 26 questions and check the judge against human grading.
- **Hosted/authenticated deployment.** Explicitly out of scope for the first
  release per [design.md](design.md#9-explicit-assumptions-mvp-scope) — no
  login or cloud dependency is planned for the near term. Ask AREE is hosted,
  but it is optional and unauthenticated, and nothing else depends on it; the
  Streamlit interface and the MCP server still run locally, and
  `aree mcp --transport streamable-http` needs authentication in front of it
  before it is shared.

## Related documentation

- [design.md](design.md) — assumptions this roadmap tracks against
- [architecture.md](architecture.md) — which directories implement which layer
