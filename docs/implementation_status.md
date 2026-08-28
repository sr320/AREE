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
| `aree` CLI | `src/aree/cli.py` | all 6 required commands run on a clean machine, plus `build-crosswalk` and `intake-supplementary` |
| Streamlit interface | `app/main.py` | browse/filter studies, evidence, candidates, search, CSV/TSV export; verified serving HTTP 200 |
| Quarto documentation site | `docs/` | 18 pages render cleanly (`quarto render docs`) |
| Reproducible supplementary intake | `src/intake/run_intake.py`, `data/studies/*/intake.yaml` | published artifact → AREE result tables via a committed config; source checksum verified before every conversion; `--check` regenerates into a temp dir and compares against committed files *and* provenance; never imputes a missing statistic |
| Species + assembly reference data | `registry/controlled_vocabularies/species.yaml`, `data/reference/genome_assemblies.yaml` | scientific names and accepted synonyms (*Crassostrea*/*Magallana gigas*) resolve to NCBI taxids; both oyster assemblies carry NCBI-verified accessions; enforced by `aree validate-study`, which also cross-checks assembly taxid against species |
| Annotation-crossing provenance | `identifier_annotation_release` in the evidence schema | records the annotation `feature_id_standardized` actually refers to, which differs from the study's own `genome_assembly` for every real study so far; empty (not invented) for the synthetic demo crosswalk |
| BioProject sample-sheet intake | `src/intake/ena_samplesheet.py`, `aree fetch-samplesheet` | reads ENA's run report + sample attributes to emit a design sheet, a FASTQ manifest with ENA MD5s, and checksummed provenance; warns when the smallest group has n<3; `validate-study` cross-checks the declared BioProject against it |
| Test suite | `tests/` (146 tests) | schema, malformed input, duplicate IDs, provenance, mapping confidence, effect sizes, meta-analysis, score reproducibility, evidence-card generation, full CLI pipeline, crosswalk construction + selection + real-crosswalk invariants, intake reproducibility + drift detection, species/assembly resolution + validation wiring |
| CI | `.github/workflows/ci.yml` | python tests + demo pipeline, real-study path (intake `--check` → register → harmonize against the real crosswalk → output assertions → clean-tree check), schema sanity, docs render |
| Lint determinism | `pyproject.toml` `[tool.ruff]` | rule set declared explicitly and ruff pinned to a minor range, so CI cannot break on an upstream default change |

## Scaffolded but not yet production-ready

| Component | Location | What's real | What's missing |
|---|---|---|---|
| RNA-seq Nextflow workflow | `workflows/rnaseq/`, `modules/rnaseq/` | **Executed end to end against real public FASTQ** (PRJNA1329250, 11 libraries, 38 processes): FASTQC → fastp → Salmon index → Salmon quant → MultiQC → DESeq2 → standardize → manifest → report. Produced 22,301 genes with real `lfcSE`. Also runs `processed_results_harmonization` end to end | run natively (Homebrew tools), **not** inside the declared containers; run on subsampled reads, not at full depth; only one of six comparisons; no methylation/proteomics/metabolomics equivalent |
| Methylation Nextflow workflow | `workflows/methylation/`, `modules/methylation/` | parses and runs `processed_results_harmonization` end to end; real Bismark/methylKit command lines | raw path never run against bisulfite data; R scripts unexecuted |
| Proteomics Nextflow workflow | `workflows/proteomics/`, `modules/proteomics/` | parses and runs `processed_results_harmonization` end to end; limma diff-abundance written | raw path not run; no raw abundance-matrix fixture |
| Metabolomics Nextflow workflow | `workflows/metabolomics/`, `modules/metabolomics/` | parses and runs `processed_results_harmonization` end to end | raw path not run against real feature-table data |
| Container images | `containers/README.md` | declared image tags per step | **no image has ever been pulled or built.** The one real execution used Homebrew-installed tools via `-profile local`, so the container specifications remain entirely unverified |

`processed_results_harmonization` is now **verified** for all four workflows —
but it was not runnable before 2026-08-28. Three of the four did not compile at
all on Nextflow 26.04 (strict DSL2 rejects script-level statements), and
metabolomics silently skipped two of its three processes while exiting 0. The
sixteen defects found and fixed are catalogued in
[first_raw_reanalysis.md](first_raw_reanalysis.md#what-the-pilot-found).

The previous version of this table claimed "DSL2 wiring correct; DAG executes
in `-stub-run` mode". That was written against an older Nextflow and had become
false; it is recorded here because a status page that quietly corrects itself is
worth less than one that says what it got wrong.

## Built but not yet exercised end-to-end

| Component | Status |
|---|---|
| Real crosswalk | Built, checksummed, unit-tested, **and exercised end-to-end on a real published study** (`HESSER2024_VCOR`): 87.2% of 351 real identifiers resolved. Now re-verified by CI on every push. |
| Meta-analysis on real evidence | The one real study contributes **no pooled estimate**: its source reports no standard error and no unadjusted p-value, so AREE declines to pool rather than impute. Random-effects pooling has therefore still only been exercised on simulated inputs. |
| Real-study coverage | Two real studies registered, both RNA-seq, both pathogen challenge. `HESSER2024_VCOR` is harmonized; `CALLA2026_OSHV` is registered but **not harmonized** — its reanalysis has not been run. Cross-study meta-analysis on real data still needs that run to happen. |
| First raw-reanalysis study | `CALLA2026_OSHV` (PRJNA1329250): design verified from deposited metadata, sample sheet + FASTQ manifest + run config complete, 226 GB of reads **not downloaded** and the workflow **not executed**. See [first_raw_reanalysis.md](first_raw_reanalysis.md). |
| Intake converter breadth | The intake config covers **differential-expression tables only** (the `gene_id`/`log2FoldChange`/`lfcSE`/`pvalue`/`padj` contract). Published methylation-region, protein-abundance, and metabolite-feature tables still need a per-study script; extending `convert_de_table` to those contracts is the obvious next increment. |

## Planned future work

- Extending the intake converter beyond differential-expression tables to the
  methylation, proteomics, and metabolomics input contracts.
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
