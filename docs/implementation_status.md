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
| 7 simulated demo studies | `registry/studies/GIGAS_*.yaml`, `registry/studies/MYTGAL_HEAT01.yaml` | RNA-seq ×2, methylation, proteomics, metabolomics, processed-only, imperfect-mapping, conflicting-direction cases all represented; a second species (`MYTGAL_HEAT01`, *M. galloprovincialis*) reuses oyster demo identifiers under the same phenotype to prove species partitioning, with oyster results tested identical with and without it. The demo list is every `simulated: true` study (`aree study-ids --simulated`) |
| Identifier harmonization | `src/harmonize/identifiers.py` | full hierarchy + 6 mapping-confidence levels incl. `unresolved`; multi-accession UniProt and multi-xref Ensembl handled; crosswalk chosen per study from its genome assembly's `identifier_crosswalk` (demo crosswalk for simulated studies), so several species harmonize in one run; `$AREE_CROSSWALK` is a taxon-checked override |
| **Real identifier crosswalks** | `data/reference/crosswalk/`, `src/mappings/build_crosswalk.py` | one per species: six mussel and clam crosswalks (41k to 49k genes each, built 2026-10-11) alongside 33,356 real *M. gigas* genes from NCBI Gene + UniProtKB, plus a 9,057-entry retired-GeneID sidecar from NCBI gene_history; checksummed provenance with per-column coverage; rebuilt by `aree build-crosswalk` |
| Per-assay harmonizers → shared evidence schema | `src/harmonize/{rnaseq,methylation,proteomics,metabolomics}.py` | both raw-reanalysis and processed-results inputs converge on one schema |
| Provenance manifests | `src/harmonize/core.py` (`_write_harmonize_manifest`) | checksum, params, versions, warnings per comparison |
| Random-effects meta-analysis | `src/meta_analysis/` | DerSimonian–Laird, I²/τ², direction-consistency; Benjamini–Hochberg `adjusted_p_value` per phenotype/feature-type/origin/species family; tail p-values via the survival function (no underflow to 0); standard errors derived once for the whole table and groups sliced as flat arrays, so a 54k-record, 30k-feature pool runs in about 5 s (was ~95 s); poolability-aware replication counts; lower-confidence alias de-duplication; fails closed on repeated within-study contrasts unless the study prespecifies one comparison as `meta_analysis_primary`, with the others counted in `n_excluded_non_primary` |
| Transparent candidate scoring + tiers | `src/prioritize/` | named weighted components scored on adjusted p, each unbounded component saturating rather than clipping so the top of a genome-wide ranking does not tie (126 candidates shared one score before 2026-09-05; the largest tie group is now 10); every candidate also carries `evidence_class` and `context_replication`, which say what the evidence is *of* and whether it replicates across biological contexts, independent of tier; hard gates that a score cannot override, including adjusted p ≤ 0.05 for the top tier; multi-omics convergence counts molecular layers (from the feature_types vocabulary) with significant same-phenotype support, not feature types anywhere; origin/species partitions preserved; tiers ordered strongest-first |
| Metabolite → gene links (`metabolite_linked_convergence` tier) | `src/mappings/metabolite_links.py`, `src/mappings/build_metabolite_links.py`, `data/mappings/metabolite_compound_map.tsv`, `data/mappings/currency_metabolites.yaml` | metabolite name → KEGG compound → reaction → KO → crg GeneID; currency and hub (> 25 genes) compounds shown but never credited; links never `exact`; separate tier below shared-identity convergence, no score contribution; demo GIGAS_GROW06 lactate ↔ simulated `ldh` converges, `sdha`/succinate and ATP do not. The real KEGG-derived table is built locally (`aree build-metabolite-links`, 191 links on 2026-10-16) and **not committed** (KEGG terms), so CI exercises the simulated table only; no real metabolomics study has passed QC yet (DUPERRET2025_POMS_STARVE failed), so no real candidate is linked |
| Evidence cards | `src/reporting/evidence_cards.py` | every card opens with a "What this evidence does and does not show" section naming the evidence class, the contexts that do and do not vary, and — where no contributing study measured a resilience outcome — the phenotype gap ahead of any validation advice; every candidate incl. emerging ranked into `candidates.tsv`; a markdown card for each with a significant signal (BH-adjusted p ≤ 0.05 pooled or in a contributing study; `--all-cards` for all); single indexed pass, so a 30k-candidate real pool renders in about a minute; partition-safe filenames, forest/multi-omics context, limitations and next step |
| `aree` CLI | `src/aree/cli.py` | all 6 required commands run on a clean machine, plus `build-crosswalk` and `intake-supplementary` |
| Streamlit interface | `app/main.py` | browse/filter studies, evidence, candidates, search, CSV/TSV export; verified serving HTTP 200 |
| Quarto documentation site | `docs/` | 18 pages render cleanly (`quarto render docs`) |
| Reproducible supplementary intake | `src/intake/run_intake.py`, `data/studies/*/intake.yaml` | published artifact → AREE result tables via a committed config; source checksum verified before every conversion; `--check` regenerates into a temp dir and compares against committed files *and* provenance; never imputes a missing statistic |
| Species + assembly reference data | `registry/controlled_vocabularies/species.yaml`, `data/reference/genome_assemblies.yaml` | scientific names and accepted synonyms (*Crassostrea*/*Magallana gigas*) resolve to NCBI taxids; both oyster assemblies and seven mussel and clam RefSeq assemblies carry NCBI-verified accessions; enforced by `aree validate-study`, which also cross-checks assembly taxid against species (a mismatch is allowed only with `reference_species_mismatch`), requires `species_identity` for species in a hybridizing complex, and warns when an assembly has no crosswalk yet |
| Annotation-crossing provenance | `identifier_annotation_release` in the evidence schema | records the annotation `feature_id_standardized` actually refers to, which differs from the study's own `genome_assembly` for every real study so far; empty (not invented) for the synthetic demo crosswalk |
| Reference bundles | `src/intake/reference_bundle.py`, `aree fetch-reference` | fetches a RefSeq transcriptome, GTF, tx2gene and optional genome for any registered assembly; every file checked against NCBI's MD5s and recorded with SHA-256 in `reference_provenance.json`, with the GTF's annotation release checked against the registry and the bundle's GeneIDs against the species' crosswalk. The RNA-seq manifest checksums the reference files a run read and compares them with the bundle. Oyster and *M. galloprovincialis* bundles fetched; the mussel one ran the full RNA-seq workflow on simulated reads (2026-10-12) |
| BioProject sample-sheet intake | `src/intake/ena_samplesheet.py`, `aree fetch-samplesheet` | reads ENA's run report + sample attributes to emit a design sheet, a FASTQ manifest with ENA MD5s, and checksummed provenance; warns when the smallest group has n<3; `validate-study` cross-checks the declared BioProject against it |
| Curated sample exclusion | `comparisons[].excluded_samples` in the study schema; `--exclude_samples` in `workflows/rnaseq` | each exclusion carries a reason and the evidence behind it; `validate-study` rejects a sample ID that is not in the sample sheet, and the workflow stops rather than excluding nothing; the run manifest records the exclusions, and a test checks that every committed CALLA2026 result was produced with exactly the exclusions its registry entry declares |
| Per-sample viral read QC | `modules/rnaseq/viral_load.nf`, `--viral_reference` | maps every read pair to a viral genome (OsHV-1 RefSeq `NC_005881.2`, committed with a provenance sidecar) and writes per-sample viral load into the manifest and report, warning on a virus-positive control or a virus-negative challenged sample. Run at full depth on two CALLA2026 comparisons, where it found an infected control (DO_CT_4, now excluded) and three uninfected challenged animals. A QC count, not a titre |
| Test suite | `tests/` (218 tests) | schema, malformed input, duplicate IDs, provenance, mapping confidence, effect sizes, poolability/repeated-contrast guards, prespecified-primary pooling, origin/species collision isolation, score reproducibility, evidence-card generation, full CLI pipeline, crosswalk construction + real-reference invariants, intake reproducibility, species/assembly validation, sample exclusions against the sample sheet and against committed manifests, committed viral-load evidence |
| CI | `.github/workflows/ci.yml` | python tests + demo pipeline, real-study path, schema sanity, docs render, all four processed-result Nextflow paths on 26.04, and a stubbed RNA-seq raw-DAG fixture test that also wires viral QC and a sample exclusion |
| Lint determinism | `pyproject.toml` `[tool.ruff]` | rule set declared explicitly and ruff pinned to a minor range, so CI cannot break on an upstream default change |

## Scaffolded but not yet production-ready

| Component | Location | What's real | What's missing |
|---|---|---|---|
| RNA-seq Nextflow workflow | `workflows/rnaseq/`, `modules/rnaseq/` | **Executed at full depth against real public FASTQ** on eight comparisons (three from PRJNA1329250, one from PRJNA593309, three from PRJNA826964 that failed QC, one from PRJNA735889): FASTQC → fastp → Salmon quant → MultiQC → DESeq2 → standardize → viral QC → manifest → report, 23,189–30,625 genes per comparison, every one with a real `lfcSE`. An arm can pool several condition levels, and `--replicate_unit` fits DESeq2 on pseudo-bulk units (e.g. tanks) where treatment was applied per unit; the manifest records both. The two 2026-09-30 CALLA2026 re-runs are the first to finish `OK`, with the report task in local scratch and the work directory off exFAT; every earlier full-depth run exited `ERR` at `RENDER_REPORT`. Processed mode runs end to end; CI also stub-runs the raw DAG | runs are native (`-profile local`, Homebrew tools), **not** inside the declared containers; manifests committed before 2026-10-11 list the declared container versions in `software_versions` (e.g. Salmon 1.10.3, DESeq2 1.42.0), not those that ran (Salmon 2.6.0, DESeq2 1.50.2, tximport 1.38.2); newer runs build it from each process's `versions.yml`, and CI checks the manifest against the work dir (`scripts/check_manifest_versions.py`); `qc_metrics` is still null; a non-stub run on the tiny CI fixtures fails in `SALMON_QUANT` under Salmon 2.6; the CI stub run validates wiring, not tools or biology |
| Methylation Nextflow workflow | `workflows/methylation/`, `modules/methylation/` | parses and runs `processed_results_harmonization` end to end; real Bismark/methylKit command lines | raw path never run against bisulfite data; R scripts unexecuted |
| Proteomics Nextflow workflow | `workflows/proteomics/`, `modules/proteomics/` | both modes run end to end in CI. Raw mode takes a wide protein x sample matrix (intensities or spectral counts/NSAF) with an explicit input scale, zeros as not detected, a per-arm detection filter and limma with posterior SE (`lfcSE`); CI runs it on a simulated spectral-count fixture (`tests/fixtures/proteomics_raw`) | no spectral search: "raw" starts from a quantified matrix, not .raw spectra |
| Metabolomics Nextflow workflow | `workflows/metabolomics/`, `modules/metabolomics/` | runs `processed_results_harmonization` end to end; CI runs the raw path for real (not stubbed) on a simulated targeted-panel fixture: per-platform internal-standard normalization, detection limit, `--replicate_unit`, limma with `lfcSE`, missingness report, manifest and report | native only (`-profile local`), no container pulled; the untargeted TIC path has no real-data run |
| Container images | `containers/README.md` | declared image tags per step | **no image has ever been pulled or built.** Every real execution used Homebrew-installed tools via `-profile local`, so the container specifications remain entirely unverified |

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
| Meta-analysis on real evidence | `CALLA2026_OSHV` and `DELISLE2020_OSHV_TEMP` pool genome-wide (23,094 shared genes, k = 2), yielding 1,994 `high_priority_cross_study` candidates. Before family-wise FDR and the significance gate were added on 2026-09-05, 12,658 of those genes qualified — 8,533 with pooled p > 0.05 — because two studies agreeing in sign was the whole test. The surviving top of the list is biologically coherent (interferon-induced helicase C domain protein, IFI44-like, rhamnose-binding lectin), which is a sanity check on the pipeline, not a validation of any candidate. `HESSER2024_VCOR` still contributes no pooled estimate: its source reports no standard error and no unadjusted p-value, so AREE declines to pool rather than impute. |
| Real-study coverage | **Updated 2026-10-09:** six real studies registered. `LUTIER2022_OA_TIPPING` (PRJNA735889, ocean acidification) is harmonized at full depth (31,425 records, fitted on 7 tanks): the first real stressor class other than pathogen challenge. It is `exposure_only` and has no OA partner, so it pools with nothing yet. `IOCAS2022_OA_ENERGY` failed QC (technical replicates) and never reaches the evidence table. `CIBNOR2021_HEAT_RRSS` is registered, with reanalysis blocked. Candidate `evidence_class` and the relevance score now follow each comparison's own `resilience_classification`, taking the most conservative across contributing records. Until 2026-10-09 they followed the phenotype default, which would have shown LUTIER's genes as `resilience_associated`. The 2026-09 summary follows. Four real studies registered. `CALLA2026_OSHV` (three comparisons, 30,625 / 30,482 / 30,541 evidence records) and `DELISLE2020_OSHV_TEMP` (23,189) are harmonized at full depth, and pool through CALLA's one prespecified comparison; `HESSER2024_VCOR` is harmonized but unpoolable; `IOCAS2022_OA_ENERGY` (PRJNA826964, ocean acidification, `acidification_tolerance`) is registered but not yet run. **The three harmonized studies are all RNA-seq pathogen challenge measuring `disease_susceptibility` or `disease_resistance` — no harmonized real study measures a resilience outcome per animal.** Every real top-tier candidate is therefore `disease_associated`, which the evidence cards state explicitly. CALLA2026's per-animal viral load is now measured from reads, but only as QC. IOCAS2022 would supply the first second stressor class once run. |
| First raw-reanalysis study | `CALLA2026_OSHV` (PRJNA1329250): three of six comparisons run at full depth and harmonized; Miyagi × USA is the `meta_analysis_primary` comparison. Viral read QC found an infected control (DO_CT_4), now excluded from every Midori comparison, and showed only 2 of 5 Midori × France animals infected. All 84 raw FASTQs (the working copies on `blue-block`) were verified against ENA's MD5s on 2026-09-30. Three comparisons remain unrun. See [first_raw_reanalysis.md](first_raw_reanalysis.md). |
| Intake converter breadth | The intake config covers **differential-expression tables only** (the `gene_id`/`log2FoldChange`/`lfcSE`/`pvalue`/`padj` contract). Published methylation-region, protein-abundance, and metabolite-feature tables still need a per-study script; extending `convert_de_table` to those contracts is the obvious next increment. |

## Planned future work

- Extending the intake converter beyond differential-expression tables to the
  methylation, proteomics, and metabolomics input contracts.
- Cross-species pooling. Orthogroups exist (OrthoFinder run
  `orthofinder3.1.5_2026-10-13`, attached to evidence at harmonization) and
  orthogroup-level meta-analysis exists (`aree meta-analyze --level
  orthogroup`), as do cross-species candidates and their evidence cards
  (`aree build-cross-species-cards`). On real data the oyster and two Manila
  clam studies pool through 5,975 single-copy orthogroups, 5,633 past the
  replication gate; none converges, because both clam comparisons are
  exposure-only. The dashboard shows cross-species results only for
  phenotypes past the WP6 gate (two species each with a pool of k >= 2;
  thermal tolerance so far) and an ortholog lookup for every top candidate;
  see WP6 in mussel_clam_integration_plan.md.
- GFF-route crosswalks for genomes without NCBI Gene records (WP7): the builder
  and *M. coruscus* (`mcoru`) exist; razor clam and geoduck need their authors'
  annotation files, and GFF-route species are not yet in an orthogroup run.
- Improving UniProt coverage beyond the 8.4% of genes that UniProtKB itself
  cross-references to NCBI GeneID — most likely by mapping RefSeq protein
  accessions through the NCBI annotation rather than relying on UniProt xrefs.
- Wiring `gene_synonyms` (already emitted in the real crosswalk) into
  `resolve_identifier` so that legacy symbols from older publications resolve.
- Additional species beyond *C. gigas* (schema already carries `species`; this
  is a data + reference addition, see [adding_a_species.md](adding_a_species.md)).
- Real-data execution of the methylation raw path. RNA-seq has tiny
  paired-FASTQ fixtures for a stubbed CI DAG test, and metabolomics and
  proteomics have fixtures CI runs for real (a targeted panel; a spectral-count
  matrix); an equivalent fixture for methylation remains future work.
- Recording observed tool versions and parsed Salmon/MultiQC QC metrics in the
  RNA-seq manifest, instead of declared versions and null `qc_metrics`.
- A covariance-aware within-study model, so a study can contribute more than
  one comparison to a pool; `meta_analysis_primary` avoids the question rather
  than answering it.
- Hosted/authenticated deployment of the interface (explicitly out of scope for
  the first release).
- Expanded phenotype/stressor ontologies and mapping to external ontology IDs.
- Genome-version liftover tooling for coordinate-based (methylation/QTL) evidence.
