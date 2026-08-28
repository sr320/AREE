# Implementation status (candid)

This is the honest Phase 4 status table required by the AREE build plan. It
distinguishes what is **complete and runnable**, what is **scaffolded but not
yet production-ready**, and what is **planned future work**. Nothing here is
sugar-coated — see [roadmap.md](roadmap.md) for the forward plan.

## Complete and runnable

Verified end-to-end on a clean install (`pytest` green, `ruff` clean, full
demo pipeline runs from the committed header-only registry).

| Component | Location | Notes |
|---|---|---|
| Study registration schema + validation | `schemas/study.schema.json`, `src/intake/` | JSON Schema + controlled-vocabulary checks; warnings vs. hard errors distinguished |
| Controlled vocabularies | `registry/controlled_vocabularies/` | phenotype (with resilience_relevance), stressor, assay, feature, tissue, life-stage, mapping-confidence, quality-flag |
| Registry ingestion | `src/intake/registry.py` | duplicate detection, `--update` path, flat CSV index |
| 6 simulated demo studies | `registry/studies/GIGAS_*.yaml` | RNA-seq ×2, methylation, proteomics, metabolomics, processed-only, imperfect-mapping, conflicting-direction cases all represented |
| Identifier harmonization | `src/harmonize/identifiers.py` | full hierarchy + 6 mapping-confidence levels incl. `unresolved`; multi-accession UniProt and multi-xref Ensembl handled; crosswalk selectable via `$AREE_CROSSWALK` |
| **Real identifier crosswalk** | `data/reference/crosswalk/`, `src/mappings/build_crosswalk.py` | 33,356 real *M. gigas* genes from NCBI Gene + UniProtKB, plus a 9,057-entry retired-GeneID sidecar from NCBI gene_history; checksummed provenance with per-column coverage; rebuilt by `aree build-crosswalk` |
| Per-assay harmonizers → shared evidence schema | `src/harmonize/{rnaseq,methylation,proteomics,metabolomics}.py` | both raw-reanalysis and processed-results inputs converge on one schema |
| Provenance manifests | `src/harmonize/core.py` (`_write_harmonize_manifest`) | checksum, params, versions, warnings per comparison |
| Random-effects meta-analysis | `src/meta_analysis/` | DerSimonian–Laird, I²/τ², direction-consistency; excludes `unresolved` from pooling |
| Transparent candidate scoring + tiers | `src/prioritize/` | named weighted components; hard gates that a score cannot override |
| Evidence cards | `src/reporting/evidence_cards.py` | one per candidate incl. emerging; forest data, multi-omics context, limitations, next step |
| `aree` CLI | `src/aree/cli.py` | all 6 required commands run on a clean machine |
| Streamlit interface | `app/main.py` | browse/filter studies, evidence, candidates, search, CSV/TSV export; verified serving HTTP 200 |
| Quarto documentation site | `docs/` | 18 pages render cleanly (`quarto render docs`) |
| Test suite | `tests/` (83 tests) | schema, malformed input, duplicate IDs, provenance, mapping confidence, effect sizes, meta-analysis, score reproducibility, evidence-card generation, full CLI pipeline, crosswalk construction + selection + real-crosswalk invariants |
| CI | `.github/workflows/ci.yml` | python tests + demo pipeline, schema sanity, docs render |

## Scaffolded but not yet production-ready

| Component | Location | What's real | What's missing |
|---|---|---|---|
| RNA-seq Nextflow workflow | `workflows/rnaseq/`, `modules/rnaseq/` | DSL2 wiring correct; DAG executes in `-stub-run` mode; real fastp/salmon/DESeq2 command lines written | not run against real FASTQ; containers not pulled; DESeq2 R script not executed |
| Methylation Nextflow workflow | `workflows/methylation/`, `modules/methylation/` | DSL2 wiring; real Bismark/methylKit command lines; both modes structured | not run against real bisulfite data; R scripts unexecuted |
| Proteomics Nextflow workflow | `workflows/proteomics/`, `modules/proteomics/` | DSL2 wiring; limma diff-abundance; Python steps stdlib/pandas-only and would run today | not run inside Nextflow; no raw abundance-matrix fixture |
| Metabolomics Nextflow workflow | `workflows/metabolomics/`, `modules/metabolomics/` | DSL2 wiring; annotation-confidence + pathway-mapping steps | not run against real feature-table data |
| Container images | `containers/README.md` | declared image tags per step | no image built or pulled in this build |

The `processed_results_harmonization` mode of each workflow depends only on
lightweight, real, pullable images (`python:3.11-slim`, Quarto) and is expected
to be runnable on a machine with Docker + Nextflow — but that has **not** been
verified here. Each workflow README carries its own "What has and has not been
run" honesty statement.

## Built but not yet exercised end-to-end

| Component | Status |
|---|---|
| Real crosswalk | Built, checksummed, unit-tested, **and now exercised end-to-end on a real published study** (`HESSER2024_VCOR`): 87.2% of 351 real identifiers resolved. |
| Meta-analysis on real evidence | The one real study contributes **no pooled estimate**: its source reports no standard error and no unadjusted p-value, so AREE declines to pool rather than impute. Random-effects pooling has therefore still only been exercised on simulated inputs. |
| Real-study coverage | One real study, one assay type (RNA-seq), one stressor (pathogen challenge). Cross-study meta-analysis on real data needs several more. |

## Planned future work

- Real ortholog mapping (OrthoFinder/OrthoDB) to populate `orthogroup_id`, which
  the real crosswalk ships empty rather than fabricating. Until then, cross-species
  evidence pooling is not supported.
- Improving UniProt coverage beyond the 8.4% of genes that UniProtKB itself
  cross-references to NCBI GeneID — most likely by mapping RefSeq protein
  accessions through the NCBI annotation rather than relying on UniProt xrefs.
- Wiring `gene_synonyms` (already emitted in the real crosswalk) into
  `resolve_identifier` so that legacy symbols from older publications resolve.
- Additional species beyond *C. gigas* (schema already carries `species`; this
  is a data + reference addition, see [adding_a_species.md](adding_a_species.md)).
- Execution of the raw-data Nextflow workflows against real public datasets,
  with small synthetic FASTQ/mzML fixtures for CI smoke-testing.
- Hosted/authenticated deployment of the interface (explicitly out of scope for
  the first release).
- Expanded phenotype/stressor ontologies and mapping to external ontology IDs.
- Genome-version liftover tooling for coordinate-based (methylation/QTL) evidence.
