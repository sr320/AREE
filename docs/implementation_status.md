# Implementation status (candid)

Reviewed against the repository and deployed data on **2026-10-18**. This page
distinguishes what is **complete and runnable**, what is **scaffolded but not
yet production-ready**, and what is **planned future work**. Nothing here is
sugar-coated — see [roadmap.md](roadmap.md) for the forward plan.

## Complete and runnable

The generated [dashboard](index.qmd) is the source for current counts.
Registry indexes are generated from study YAMLs, not committed. CI exercises
the demo pipeline and the HESSER2024 published-source intake path.

| Component | Location | Notes |
|---|---|---|
| Study registration schema + validation | `schemas/study.schema.json`, `src/intake/` | JSON Schema + controlled-vocabulary checks; warnings vs. hard errors distinguished |
| Controlled vocabularies | `registry/controlled_vocabularies/` | phenotype (with resilience_relevance), stressor, assay, feature, tissue, life-stage, mapping-confidence, quality-flag |
| Registry ingestion | `src/intake/registry.py` | duplicate detection, `--update` path, flat CSV index |
| 7 simulated demo studies | `registry/demo_studies/GIGAS_*.yaml`, `registry/demo_studies/MYTGAL_HEAT01.yaml` | RNA-seq ×4 (including the mussel and legacy-ID study), methylation, proteomics, metabolomics, processed-only, imperfect-mapping, conflicting-direction cases all represented; a second species (`MYTGAL_HEAT01`, *M. galloprovincialis*) reuses oyster demo identifiers under the same phenotype to prove species partitioning, with oyster results tested identical with and without it. The demo list is every `simulated: true` study (`aree study-ids --simulated`) |
| Identifier harmonization | `src/harmonize/identifiers.py` | full hierarchy + 6 mapping-confidence levels incl. `unresolved`; multi-accession UniProt and multi-xref Ensembl handled; crosswalk chosen per study from its genome assembly's `identifier_crosswalk` (demo crosswalk for simulated studies), so several species harmonize in one run; `$AREE_CROSSWALK` is a taxon-checked override |
| **Real identifier crosswalks** | `data/reference/crosswalk/`, `src/mappings/build_crosswalk.py` | one per species: six mussel and clam crosswalks (41k to 49k genes each, built 2026-10-11) alongside 33,356 real *M. gigas* genes from NCBI Gene + UniProtKB, plus a 9,057-entry retired-GeneID sidecar from NCBI gene_history; checksummed provenance with per-column coverage; rebuilt by `aree build-crosswalk` |
| Per-assay harmonizers → shared evidence schema | `src/harmonize/{rnaseq,methylation,proteomics,metabolomics}.py` | both raw-reanalysis and processed-results inputs converge on one schema |
| Provenance manifests | `src/harmonize/core.py` (`_write_harmonize_manifest`) | checksum, params, versions, warnings per comparison |
| Random-effects meta-analysis | `src/meta_analysis/` | DerSimonian–Laird, I²/τ², direction-consistency; Benjamini–Hochberg `adjusted_p_value` per phenotype/feature-type/origin/species/contrast-type family; tail p-values via the survival function (avoids cancellation; extreme tails can still underflow); standard errors derived once for the whole table and groups sliced as flat arrays, so a 54k-record, 30k-feature pool runs in about 5 s (was ~95 s); poolability-aware replication counts; lower-confidence alias de-duplication; fails closed on repeated within-study contrasts unless the study prespecifies one comparison as `meta_analysis_primary`, with the others counted in `n_excluded_non_primary` |
| Transparent candidate scoring + tiers | `src/prioritize/` | named weighted components scored on adjusted p, each unbounded component saturating rather than clipping so the top of a genome-wide ranking does not tie (126 candidates shared one score before 2026-09-05 in the initial pool); every candidate also carries `evidence_class` and `context_replication`, which say what the evidence is *of* and whether it replicates across biological contexts, independent of tier; hard gates that a score cannot override, including adjusted p ≤ 0.05 for the top tier; multi-omics convergence counts molecular layers (from the feature_types vocabulary) with significant same-phenotype support, not feature types anywhere; origin/species partitions preserved; tiers ordered strongest-first |
| Metabolite → gene links (`metabolite_linked_convergence` tier) | `src/mappings/metabolite_links.py`, `src/mappings/build_metabolite_links.py`, `data/mappings/metabolite_compound_map.tsv`, `data/mappings/transporter_substrates.tsv`, `data/mappings/currency_metabolites.yaml` | metabolite name → KEGG compound → reaction → KO → crg GeneID (enzymes), plus curated transporter KOs verified against KEGG BRITE ko02000 (transporters, always `inferred`; 191 enzyme + 58 transporter links on 2026-10-17); currency and hub (> 25 genes) compounds shown but never credited; links never `exact`; separate tier below shared-identity convergence, no score contribution; demo GIGAS_GROW06 lactate ↔ simulated `ldh` (enzyme) and glucose ↔ simulated `glut1` (transporter) converge, `sdha`/succinate and ATP do not. The real KEGG-derived table is built locally (`aree build-metabolite-links`) and **not committed** (KEGG terms), so CI exercises the simulated table only; no real metabolomics study has passed QC yet (DUPERRET2025_POMS_STARVE failed), so no real candidate is linked |
| Evidence cards | `src/reporting/evidence_cards.py` | every card opens with a "What this evidence does and does not show" section naming the evidence class, the contexts that do and do not vary, and — where no contributing study measured a resilience outcome — the phenotype gap ahead of any validation advice; every candidate incl. emerging ranked into `candidates.tsv`; a markdown card for each with a significant signal (BH-adjusted p ≤ 0.05 pooled or in a contributing study; `--all-cards` for all); single indexed pass, so a 30k-candidate real pool renders in about a minute; partition-safe filenames, forest/multi-omics context, limitations and next step |
| `aree` CLI | `src/aree/cli.py` | all 6 required commands run on a clean machine, plus `build-crosswalk` and `intake-supplementary` |
| Streamlit interface | `app/main.py` | browse/filter studies, evidence, candidates, search, CSV/TSV export; verified serving HTTP 200 |
| Quarto documentation site | `docs/` | 32 documentation source pages, plus generated evidence-card pages (`quarto render docs`) |
| Reproducible supplementary intake | `src/intake/run_intake.py`, `data/studies/*/intake.yaml` | published artifact → AREE result tables via a committed config; source checksum verified before every conversion; `--check` regenerates into a temp dir and compares against committed files *and* provenance; never imputes a missing statistic |
| Per-species views (WP8) | `src/reporting/species.py`, `docs/index.qmd`, `app/main.py` | one row per registered species with its crosswalk route and size, real and simulated study counts and evidence; the dashboard's species table and studies filter, and the app's species filters on studies, evidence and candidates, all on accepted names so synonyms group; a species-identity column for studies of a hybridizing complex; `species_identity_basis` in the registry index ([species_identity.md](species_identity.md)) |
| Species + assembly reference data | `registry/controlled_vocabularies/species.yaml`, `data/reference/genome_assemblies.yaml` | scientific names and accepted synonyms (*Crassostrea*/*Magallana gigas*) resolve to NCBI taxids; both oyster assemblies and seven mussel and clam RefSeq assemblies carry NCBI-verified accessions; enforced by `aree validate-study`, which also cross-checks assembly taxid against species (a mismatch is allowed only with `reference_species_mismatch`), requires `species_identity` for species in a hybridizing complex, and warns when an assembly has no crosswalk yet |
| Annotation-crossing provenance | `identifier_annotation_release` in the evidence schema | records the annotation `feature_id_standardized` actually refers to, which may differ from the study's own `genome_assembly`; empty (not invented) for the synthetic demo crosswalk |
| Reference bundles | `src/intake/reference_bundle.py`, `aree fetch-reference` | fetches a RefSeq transcriptome, GTF, tx2gene and optional genome for any registered assembly; every file checked against NCBI's MD5s and recorded with SHA-256 in `reference_provenance.json`, with the GTF's annotation release checked against the registry and the bundle's GeneIDs against the species' crosswalk. The RNA-seq manifest checksums the reference files a run read and compares them with the bundle. Oyster and *M. galloprovincialis* bundles fetched; the mussel one ran the full RNA-seq workflow on simulated reads (2026-10-12) |
| BioProject sample-sheet intake | `src/intake/ena_samplesheet.py`, `aree fetch-samplesheet` | reads ENA's run report + sample attributes to emit a design sheet, a FASTQ manifest with ENA MD5s, and checksummed provenance; warns when the smallest group has n<3; `validate-study` cross-checks the declared BioProject against it |
| Curated sample exclusion | `comparisons[].excluded_samples` in the study schema; `--exclude_samples` in `workflows/rnaseq` | each exclusion carries a reason and the evidence behind it; `validate-study` rejects a sample ID that is not in the sample sheet, and the workflow stops rather than excluding nothing; the run manifest records the exclusions, and a test checks that every committed CALLA2026 result was produced with exactly the exclusions its registry entry declares |
| Per-sample viral read QC | `modules/rnaseq/viral_load.nf`, `--viral_reference` | maps every read pair to a viral genome (OsHV-1 RefSeq `NC_005881.2`, committed with a provenance sidecar) and writes per-sample viral load into the manifest and report, warning on a virus-positive control or a virus-negative challenged sample. Run at full depth on two CALLA2026 comparisons, where it found an infected control (DO_CT_4, now excluded) and three uninfected challenged animals. A QC count, not a titre |
| Test suite | `tests/` (use `pytest --collect-only -q` for the current case count) | schema, malformed input, duplicate IDs, provenance, mapping confidence, effect sizes, poolability/repeated-contrast guards, prespecified-primary pooling, origin/species collision isolation, score reproducibility, evidence-card generation, full CLI pipeline, crosswalk construction + real-reference invariants, intake reproducibility, species/assembly validation, sample exclusions against the sample sheet and against committed manifests, committed viral-load evidence |
| CI | `.github/workflows/ci.yml` | python tests + demo pipeline, real-study path, schema sanity, docs render, all four processed-result Nextflow paths on 26.04, and a stubbed RNA-seq raw-DAG fixture test that also wires viral QC and a sample exclusion |
| Lint determinism | `pyproject.toml` `[tool.ruff]` | rule set declared explicitly and ruff pinned to a minor range, so CI cannot break on an upstream default change |

## Scaffolded but not yet production-ready

| Component | Location | What's real | What's missing |
|---|---|---|---|
| RNA-seq Nextflow workflow | `workflows/rnaseq/`, `modules/rnaseq/` | Full-depth real public FASTQ runs across multiple oyster and Manila clam studies; FastQC → fastp → Salmon → MultiQC → DESeq2 → standardize → optional viral QC → manifest → report. Supports pooled condition levels, replicate units, curated exclusions and 3′ Tag-seq. See the registry and committed manifests for each run | native execution; declared containers unverified. Manifests before 2026-10-11 record declared rather than observed versions; newer runs collect process `versions.yml`. Parsed general QC metrics remain incomplete; CI stub runs verify wiring rather than biology |
| Methylation Nextflow workflow | `workflows/methylation/`, `modules/methylation/` | Real `CIBNOR2021_HEAT_RRSS` bisulfite run completed on 2026-10-16: Bismark plus family-level pooling and limma, 132,587 harmonized regions, no BH-significant regions. Processed mode and configurable methylKit/limma paths exist | this run supports the Bismark + pooled-limma path, not every methylKit mode or EM-seq design; declared containers unverified |
| Proteomics Nextflow workflow | `workflows/proteomics/`, `modules/proteomics/` | Both modes exercised in CI. Quantified-matrix mode uses detection filtering and limma with posterior SE. Real `TIMMINS2014_OA_GILL` results contribute 1,896 evidence records across three comparisons | no spectral search: “raw” starts from a quantified matrix, not instrument spectra; no significant records in this study’s current results |
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

## Real-data coverage and limitations

| Component | Status |
|---|---|
| Real crosswalk | Built, checksummed, unit-tested, **and exercised end-to-end on a real published study** (`HESSER2024_VCOR`): 87.2% of 351 real identifiers resolved. Now re-verified by CI on every push. |
| Meta-analysis on real evidence | Cross-study gene pools exist for oyster disease susceptibility (up to five studies), oyster thermal tolerance (up to three), and Manila clam thermal tolerance (two). Baseline and exposure-response contrasts remain separate. Current counts and evidence classes are on the [dashboard](index.qmd). `HESSER2024_VCOR` remains unpoolable because its source has neither SE nor unadjusted p-value; it remains visible as unscored emerging evidence. |
| Real-study coverage | Snapshot on 2026-10-18: **30 real studies**, 69 registered comparisons; **16 studies / 40 comparisons** with harmonized evidence, totaling **1,066,440 records** across RNA-seq, methylation and proteomics in oyster and Manila clam. Five studies have failed QC and are excluded from the Pages harmonization run. No real metabolomics study has passed QC. These are dated counts; use the generated dashboard for later builds. |
| First raw-reanalysis study | `CALLA2026_OSHV` (PRJNA1329250): three of six comparisons run at full depth and harmonized; Miyagi × USA is the `meta_analysis_primary` comparison. Viral read QC found an infected control (DO_CT_4), now excluded from every Midori comparison, and showed only 2 of 5 Midori × France animals infected. All 84 raw FASTQs (the working copies on `blue-block`) were verified against ENA's MD5s on 2026-09-30. Three comparisons remain unrun. See [first_raw_reanalysis.md](first_raw_reanalysis.md). |
| Intake converter breadth | The intake config covers **differential-expression tables only** (the `gene_id`/`log2FoldChange`/`lfcSE`/`pvalue`/`padj` contract). Published methylation-region, protein-abundance, and metabolite-feature tables still need a per-study script; extending `convert_de_table` to those contracts is the obvious next increment. |

## Planned future work

- Extending the intake converter beyond differential-expression tables to the
  methylation, proteomics, and metabolomics input contracts.
- Broader cross-species replication. Orthogroups exist (OrthoFinder run
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
- Real mussel studies (none harmonized yet; the demo mussel is simulated), and
  a clam study with a measured tolerance outcome in its pooled comparison.
- Independent methylation replication and additional raw-path validation,
  including the registered EM-seq study. CIBNOR2021 has completed a real
  Bismark/pooled-limma run; this does not validate every methylKit configuration.
- Parsing Salmon/MultiQC QC metrics into RNA-seq manifests. Observed tool
  versions are already collected for newer runs; old manifests retain the
  declared-version limitation described above.
- A covariance-aware within-study model, so a study can contribute more than
  one comparison to a pool; `meta_analysis_primary` avoids the question rather
  than answering it.
- Hosted/authenticated deployment of the interface (explicitly out of scope for
  the first release).
- Expanded phenotype/stressor ontologies and mapping to external ontology IDs.
- Genome-version liftover tooling for coordinate-based (methylation/QTL) evidence.
