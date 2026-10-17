# Plan: integrating mussel and clam datasets

**Status, 2026-10-17:** WP0 to WP4 are implemented. WP5 has its backlog and
its first within-species pool: two Manila clam heat studies
(TUCCI2025_HEAT_PRIMING and PERUZZA2023_HEATWAVE, both reanalyzed and on
`main`, k = 2). WP6 has its
orthogroup run, orthogroup-level meta-analysis and cross-species candidates
with their evidence cards, and the dashboard gate (see the entries in section 5). Everything else is a proposal. This plan describes the work needed
to take AREE from a single-species (*Magallana
gigas*) evidence engine to one that holds mussel and clam evidence, first
species-by-species and then across species through orthology.

The reference-genome and public-data facts below were queried from the NCBI
Datasets API and the ENA portal API on 2026-10-11. Re-verify them before
acting on any one of them, as [candidate_studies.md](candidate_studies.md)
requires for every accession.

## 1. What "integration" means here

Three distinct capabilities, in order of difficulty:

1. **Within-species evidence for a new species.** Register a *Mytilus* or
   *Ruditapes* study, reanalyze or harmonize it, and get species-specific
   meta-analysis, candidates and evidence cards. The existing partition on
   `species_taxid` means this never pools with oyster evidence.
2. **Several species in one run.** Harmonize oyster, mussel and clam studies
   in the same `make real-pool` / Pages build without the crosswalk, the
   reference bundle, or the dashboard assuming oysters.
3. **Cross-species convergence.** Ask whether the same gene family is
   associated with thermal tolerance in oysters and mussels. This needs an
   orthology layer that does not exist yet, and a deliberate decision about
   when pooling across species is scientifically defensible at all.

Capabilities 1 and 2 are mostly plumbing and are the MVP of this plan.
Capability 3 is new science and is scoped as its own work package with its
own gates.

## 2. What already works, and what is hard-wired to oysters

[adding_a_species.md](adding_a_species.md) is accurate: `species` is a data
field, meta-analysis and ranking partition on `species_taxid`, and
`aree validate-study` cross-checks species against the assembly's taxid. The
crosswalk builder is already parameterized by taxid and output slug
(`aree build-crosswalk --taxid --slug --organism`).

What was not species-agnostic before WP1, found by reading the code rather than
the docs. WP0 and WP1 have since resolved every row of this table except the
reference bundles, which are WP3:

| Where | What is hard-wired | Why it blocks a mixed-species pool |
|---|---|---|
| `src/harmonize/identifiers.py`, `src/harmonize/core.py:79` | One active crosswalk for the whole process, chosen by `simulated` plus `$AREE_CROSSWALK`. `crosswalk_for_study()` takes only a boolean. | A mussel study harmonized under the oyster crosswalk resolves nothing (or worse, resolves a colliding `LOC` number to an oyster gene). |
| `Makefile:64`, `.github/workflows/pages.yml:63`, `.github/workflows/ci.yml:107` | `AREE_CROSSWALK=data/reference/crosswalk/mgigas_gene_id_crosswalk.tsv` | `make real-pool` and the Pages build harmonize every committed real result against the oyster crosswalk. |
| `src/reporting/dashboard.py:36` | `REAL_CROSSWALK` points at the mgigas file | Dashboard crosswalk statistics would describe only oysters. |
| `scripts/check_real_study_evidence.py:76` | Asserts every real record has `species_taxid == 29159` | Fails the moment a real mussel record exists. This is a guard that did its job; it needs to become a per-study check. |
| `data/reference/<ASSEMBLY>/` layout, `config/<STUDY>.config` | `transcript_fasta`, `gtf`, `tx2gene` paths per study; the bundle recipe in [first_raw_reanalysis.md](first_raw_reanalysis.md) is written for one NCBI RefSeq assembly | Works per study already, but there is no command or checklist for building a second species' bundle, and the methylation workflow needs a genome FASTA plus Bismark index that the oyster bundle does not include. |
| `registry/controlled_vocabularies/tissue_types.yaml` | Nine bivalve terms, none for foot, byssal gland, siphon or hemocytes | Mussel byssus and clam siphon studies would be forced onto `mixed_tissue`. |
| `stressor_ontology.yaml`, `phenotype_ontology.yaml` | No emersion/desiccation stressor; no byssal attachment or burrowing phenotype | Intertidal mussel and infaunal clam resilience traits have no term. |
| `data/mappings/gene_id_crosswalk.tsv` (demo), `tests/` | Single-species synthetic crosswalk; species-partition tests use oyster vs. a fabricated collision | No test exercises two real-shaped species in one evidence table. |

## 3. Which species, and why

The decisive constraint is **where gene identifiers come from**. AREE's
crosswalk is built from NCBI Gene records for a taxid, and NCBI Gene only has
records for genomes that NCBI has annotated (RefSeq, `GCF_` accessions). For a
species with only a GenBank (`GCA_`) assembly and an author-supplied GFF, the
builder has nothing to download, and the whole identifier hierarchy in
[identifier_mapping.md](identifier_mapping.md) has to be rebuilt from the GFF
plus UniProt. That is a separate, larger piece of work (section 5, WP7).

### 3.1 Reference genome and public data survey (2026-10-11)

NCBI RefSeq-annotated genomes exist for nine mussel and clam species. Gene
counts are from the NCBI annotation report. ENA counts are `read_run` records
under `tax_tree(taxid)` by `library_strategy`, and "studies" is the number of
distinct study accessions. Counts include mislabelled and non-resilience
studies, so treat them as scale, not a backlog.

| Species | NCBI taxid | RefSeq assembly (annotation) | Genes (protein-coding) | ENA RNA-Seq runs / studies | Bisulfite runs / studies | WGS runs / studies |
|---|---|---|---|---|---|---|
| *Mytilus galloprovincialis* | 29158 | GCF_965363235.1 xbMytGall1.hap1.1, chromosome (RS_2025_08) | 43,856 (27,168) | 716 / 47 | 0 | 337 / 15 |
| *Mytilus edulis* | 6550 | GCF_963676685.1 xbMytEdul2.2, chromosome (RS_2025_02) | 49,209 (35,875) | 783 / 39 | 0 | 166 / 7 |
| *Mytilus trossulus* | 6551 | GCF_036588685.1 PNRI_Mtr1.1.1.hap1, chromosome (RS_2024_02) | 47,780 (34,235) | 54 / 8 | 0 | 2 / 2 |
| *Mytilus californianus* | 6549 | GCF_021869535.1 xbMytCali1.0.p, scaffold (AR 100, 2022-11) | 41,296 (26,243) | 168 / 10 | **94 / 2** | 19 / 4 |
| *Perna viridis* | 73031 | GCF_037379345.1 ASM3737934v1, chromosome (RS_2026_09) | 33,914 (24,439) | 98 / 14 | 0 | 8 / 2 |
| *Mercenaria mercenaria* | 6596 | GCF_021730395.1 MADL_Memer_1, chromosome (RS_2023_01); a second RefSeq-annotated assembly GCF_014805675.1 also exists | 45,375 (38,330) | 263 / 23 | 0 | 133 / 3 |
| *Ruditapes philippinarum* | 129788 | GCF_026571515.1 ASM2657151v2, **contig-level, N50 183 kb** (RS_2023_10); several unannotated chromosome-level GCA assemblies | 48,037 (40,909) | **1,085 / 75** | 0 | 393 / 14 |
| *Mya arenaria* | 6604 | GCF_026914265.1 ASM2691426v1, chromosome (RS_2023_01) | 44,373 (36,031) | 38 / 4 | 2 / 1 | 22 / 1 |
| *Dreissena polymorpha* (freshwater, not aquaculture) | 45954 | GCF_020536995.1 | | 109 | 0 | 177 |

Species with public data but **no RefSeq annotation** (GenBank assembly only,
or none):

| Species | NCBI taxid | Best assembly | ENA RNA-Seq runs / studies | Note |
|---|---|---|---|---|
| *Mytilus coruscus* | 42192 | GCA_059745995.1 chromosome; GenBank annotation on older MCOR1.1 | 383 / 40 | Major Chinese aquaculture mussel; large literature on hypoxia, heat, and selective breeding |
| *Sinonovacula constricta* (razor clam) | 98310 | GCA_007844125.1 chromosome, no NCBI annotation | 416 / 26 | Chinese aquaculture; breeding lines for hypoxia and ammonia tolerance |
| *Tegillarca granosa* (blood clam) | 220873 | GCA_029721355.1 chromosome, GenBank annotation | 225 / 24 | |
| *Panopea generosa* (geoduck) | 1049056 | GCA_029582155.1 chromosome; PGEN-v1.0 scaffold, no NCBI annotation | 159 / 6 | Roberts Lab species; OA and epigenetics work. No ENA runs labelled Bisulfite-Seq were found under this taxid, so the lab's methylation deposits need a label check |
| *Meretrix meretrix* | 291251 | GCA_978046135.1 chromosome | 91 / 6 | |
| *Ruditapes decussatus* | 104385 | GCA_977926895.1 chromosome | 74 / 5 | European grooved carpet shell; Perkinsus and brown ring disease |
| *Mytilus chilensis* | 173667 | **none** | 38 / 4 | Reads would have to map to *M. galloprovincialis* or *M. edulis* |
| *Perna canaliculus* (greenshell) | 38949 | **none** | 21 / 3 | NZ aquaculture; summer mortality; processed-results mode only |

### 3.2 Recommendation

**First wave (RefSeq route, no new identifier machinery):**

- Mussels: *M. galloprovincialis*, *M. edulis*, *M. trossulus*. These three
  are the *M. edulis* species complex. Together they have the most public
  RNA-seq, chromosome-level Sanger/PNRI assemblies with 2024 to 2025 RefSeq
  annotations, and a resilience literature (thermal tolerance, mass
  mortality, Vibrio challenge) that maps onto existing AREE phenotype terms.
  *M. californianus* follows because it carries the only sizeable mussel
  bisulfite deposits.
- Clams: *R. philippinarum* (Manila clam) and *M. mercenaria* (hard clam).
  Manila clam is the single largest public RNA-seq resource in the survey and
  has resistant-versus-susceptible designs for brown ring disease and
  *Perkinsus*. Hard clam has QPX-resistance lines and multi-stressor
  warming/hypoxia/acidification work, and a chromosome-level RefSeq genome.

**Second wave (needs WP7, GFF-based crosswalks):** *M. coruscus*,
*S. constricta*, *P. generosa*. Geoduck is a lab priority, but putting it
first would mean building the non-RefSeq identifier path before any
mussel or clam evidence exists. Until WP7 lands, geoduck and razor clam
studies can still be **registered** and their **processed results
harmonized with `unresolved` identifiers**, which keeps the registry and
study-quality views honest without pretending the genes are mapped.

**Deferred:** *M. chilensis*, *P. canaliculus* (no genome), *Dreissena*
(not aquaculture).

## 4. Scientific issues specific to mussels and clams

These shape the data model changes in section 5 and should be settled before
the first mussel study is registered.

**Species identity in the *Mytilus edulis* complex.** *M. edulis*,
*M. galloprovincialis* and *M. trossulus* hybridize across their contact
zones, and public deposits are labelled by the collecting lab's assumption,
not by genotyping. A study labelled *M. edulis* from a Scottish or Dutch site
may be hybrid, and a "*M. galloprovincialis*" deposit from Puget Sound is
probably *M. trossulus*. AREE must not silently trust the label. Implemented in WP0:
a `species_identity` block in the study YAML (`basis: genotyped | morphology |
collection_site | authors_statement | unstated`, plus notes), the quality flag
`species_assignment_unverified`, and a `complex: mytilus_edulis_complex` field
on the three vocabulary terms. `aree validate-study` requires the block for
complex members and requires the flag whenever the basis is not `genotyped`.
The flag is declared by the curator and enforced, rather than added silently,
so the YAML always says what the evidence records will carry.
Partitioning still uses `species_taxid`; the complex field is for display and
for WP6's orthogroup view, never for pooling by default.

**Reference-species mismatch.** A *M. chilensis* study mapped to the
*M. galloprovincialis* genome, or a geoduck study mapped to nothing, is a
different mapping situation from an oyster study on an old assembly. Add a
quality flag `reference_species_mismatch` when the study's species taxid
differs from the assembly's taxid, and require `aree validate-study` to accept
that pair only when the flag is present. Implemented in WP0; without the flag
the pair still fails outright, which is the right default.

**Contig-level RefSeq for Manila clam.** The RefSeq-annotated assembly has an
N50 of 183 kb while several unannotated chromosome-level assemblies exist.
Gene-level evidence (RNA-seq, proteomics) is unaffected because it keys on
NCBI GeneID. Coordinate-based evidence (methylation regions) on this assembly
will be fragmented and will not lift over cleanly to a future chromosome-level
RefSeq. Record this in the assembly's `notes` and apply the existing
genome-version caveats from [handling_genome_versions.md](handling_genome_versions.md).

**Two RefSeq-annotated hard clam assemblies.** NCBI has annotated both
GCF_021730395.1 (current reference) and GCF_014805675.1. The crosswalk is
built per taxid, not per assembly, so GeneIDs are shared; the retired-GeneID
sidecar handles IDs that changed between annotation releases. Register both
assemblies in `genome_assemblies.yaml` with only one marked
`is_ncbi_reference`, which the existing test enforces.

**New phenotypes and stressors.** Mussel resilience traits that oysters do
not have: byssal attachment strength and thread production (a direct fitness
and farm-retention trait under heat and acidification), and condition under
emersion. Clam traits: burrowing speed and depth after stress, and siphon
regeneration. Proposed additions, each with `resilience_relevance`:

- phenotype `byssal_attachment` (resilience), `burrowing_performance`
  (resilience), `shell_repair` (recovery, optional).
- stressor `emersion_desiccation` (air exposure at low tide or during
  transport, distinct from `handling_transport_stress`).
- tissue `foot`, `byssal_gland`, `siphon`, `hemocytes` (hemocytes are the
  cell fraction of `hemolymph` and are what most immune RNA-seq actually
  profiles; keep both terms rather than overloading one).

Keep the list to terms a first-wave study actually needs; add others when a
study needs them, as the vocabularies' own guidance says.

**Transmissible neoplasia.** Bivalve transmissible neoplasia (BTN) occurs in
*M. trossulus*, *M. edulis*, *M. chilensis* and *Mya arenaria*, and part of
the public RNA-seq for these species is cancer biology, not resilience. It
maps onto `disease_susceptibility` with `pathogen_challenge` as the stressor
only loosely. Register such studies with `resilience_classification:
disease_associated` and say in `limitations` that the "pathogen" is a
transmissible cell line, so they are not pooled with Vibrio challenges by
accident. Do not add a BTN-specific term until there are two such studies.

## 5. Work packages

Each work package is sized for a single PR on its own branch, per
[parallel_study_runs.md](parallel_study_runs.md): code and vocabulary changes
merge first, study branches rebase onto them.

### WP0. Vocabularies and reference data

**Status: done, 2026-10-11.** Six species terms (the first wave is five
species, plus *M. californianus*) and seven assemblies, including both
RefSeq-annotated hard clam assemblies. Every accession and synonym was
re-verified against NCBI on that date. WP0 needed some code after all: the
species-complex and reference-mismatch validation rules, the
`species_identity` schema field, and tests in `tests/test_multi_species.py`.

Files: `registry/controlled_vocabularies/species.yaml`,
`tissue_types.yaml`, `phenotype_ontology.yaml`, `stressor_ontology.yaml`,
`data/reference/genome_assemblies.yaml`, `tests/test_reference_data.py`.

- Add species terms for the five first-wave species and *M. californianus*,
  each with `ncbi_taxid`, the synonyms NCBI Taxonomy lists for that taxid
  (*Venerupis*, *Tapes* and *Venus philippinarum* for Manila clam; *Venus
  mercenaria* for hard clam; none for the mussels), and the `complex` field
  for the three *Mytilus*.
- Add one `genome_assemblies.yaml` entry per RefSeq assembly in section 3.1,
  with the verification date, plus the second hard clam assembly.
- Add the tissue, phenotype and stressor terms from section 4.
- Extend the reference-data tests so every species term has exactly one
  reference assembly and every assembly's taxid resolves to a species term.

Acceptance: `aree validate-study` passes on a stub *M. galloprovincialis*
YAML and on a stub Manila clam YAML; all existing oyster studies still pass.

### WP1. Per-study crosswalk selection

**Status: done, 2026-10-11.** Each assembly carries `identifier_crosswalk`, and
a real study resolves against its assembly's. `$AREE_CROSSWALK` is honoured
only for a crosswalk whose provenance taxid matches. Implementing this exposed
a latent bug: the crosswalk chosen for one study stuck for every later study in
the same process, so a mixed-species pool would have resolved mussel IDs
against the oyster crosswalk. A regression test now covers it. Harmonizing all
committed real oyster studies gives an evidence table identical to the
previous code's, apart from the absolute checkout path in `source_file`.

Three guards planned for WP2 were pulled forward, because the updated
`adding_a_species.md` now advertises `aree build-crosswalk --taxid`. The
command refuses a non-oyster taxid without `--slug`, which would have
overwritten the oyster crosswalk. It refuses to overwrite a crosswalk built for
another taxid. It derives the provenance organism and annotation release from
the registered reference assembly instead of defaulting to the oyster's.

Files: `src/harmonize/identifiers.py`, `src/harmonize/core.py`,
`src/common/__init__.py`, `Makefile`, `.github/workflows/{ci,pages}.yml`,
`src/reporting/dashboard.py`, `scripts/check_real_study_evidence.py`,
`scripts/harmonize_committed_results.py`, `tests/`.

- Add a `crosswalk` field to each `genome_assemblies.yaml` entry (or a
  taxid-keyed `data/reference/crosswalk/index.yaml`; the assembly file is
  simpler and already validated). The oyster entries point at
  `mgigas_gene_id_crosswalk.tsv`.
- Change `crosswalk_for_study(simulated: bool)` to
  `crosswalk_for_study(study: dict)`: simulated studies use the demo
  crosswalk; real studies resolve species, then assembly, then the registered
  crosswalk path. A real study whose species has no committed crosswalk is an
  error that names the `aree build-crosswalk` invocation to fix it, as today.
- `$AREE_CROSSWALK` becomes an explicit override for local experiments and is
  removed from the Makefile and both workflows. `make real-pool` and the
  Pages build then need no per-species edits, which is the same property
  PR #13 established for the study list.
- Dashboard crosswalk statistics become per-species rows.
- `check_real_study_evidence.py` checks each record's `species_taxid` against
  its own study's species rather than against 29159.
- Tests: a two-species evidence table built from fixtures with colliding
  `LOC` numbers must resolve each to its own crosswalk and produce two
  partitions with no shared `feature_id_standardized`; harmonizing a real
  mussel stub without a mussel crosswalk must fail loudly.

Acceptance: `make real-pool` with one oyster and one mussel study committed
harmonizes both, and the evidence table's `identifier_crosswalk` provenance
names a different file and checksum for each.

### WP2. Crosswalks for the first-wave species

**Status: done, 2026-10-11.** Six crosswalks built and registered, one per
species; coverage is in [identifier_mapping.md](identifier_mapping.md#mussel-and-clam-crosswalks).
Findings that change later work:

- **UniProt linkage is near zero** for every mussel and clam (2 to 94 genes per
  species), because UniProtKB entries for them derive from GenBank, not RefSeq.
  Proteomics studies (WP5) will need sequence-based protein mapping before they
  can contribute.
- **Hard clam:** the two annotated assemblies share 31,226 GeneIDs. On the older
  one, 3,656 GeneIDs remap to a current gene and 9,078 were discontinued without
  replacement, so both assemblies use one crosswalk and the older one carries a
  note that up to about 21% of its genes will stay unresolved.
- NCBI's `gene_info` and `gene_history` were confirmed to keep each taxid's rows
  contiguous, which the builder's early-exit streaming depends on.
- The real-crosswalk invariant tests now run over every registered crosswalk,
  and a new test checks that no GeneID appears in two species' crosswalks.

Run `aree build-crosswalk --taxid <id> --slug <slug>` for each of the six
species (five first wave plus *M. californianus*), commit the TSV and
provenance sidecar as for oysters (3 MB each), and set `identifier_crosswalk`
on each species' assemblies. Provenance organism and annotation release now
default to the registered reference assembly (done in WP1). For hard clam,
confirm whether NCBI Gene carries one GeneID set across both annotated
assemblies before pointing both at one crosswalk. Extend
`tests/test_crosswalk.py` real-reference invariants to loop over every
committed crosswalk (one row per GeneID, empty orthogroups, honest coverage).

Acceptance: six new crosswalks committed with provenance; test suite green.

### WP3. Reference bundles and workflow configs

**Status: done, 2026-10-12.** `aree fetch-reference` fetches, verifies and
records a bundle for any registered RefSeq assembly; see
[reference_bundles.md](reference_bundles.md). What was done and found:

- **Every registered RefSeq assembly resolves** to its NCBI folder (`--plan`),
  all eight of them; the Roslin oyster assembly, which has no RefSeq
  annotation, is refused.
- **The oyster bundle was backfilled** with `reference_provenance.json` and the
  genome FASTA that `config/CIBNOR2021_HEAT_RRSS.config` already expected. The
  new tx2gene builder reproduces the pilot's awk-built `tx2gene.tsv` byte for
  byte (60,992 transcripts, 33,068 genes, all in the oyster crosswalk).
- **The *M. galloprovincialis* bundle was fetched.** Its GTF names the
  registered release (`GCF_965363235.1-RS_2025_08`), and all 43,185 of its
  genes are in the WP2 crosswalk.
- **Acceptance, run with simulated reads.** No mussel study is curated yet
  (WP5), so six single-end libraries (3 vs 3, about 85,000 to 97,000 reads
  each) were simulated from 400 single-isoform mussel genes, 60 of them with a
  4-fold change. The real RNA-seq workflow ran end to end on this machine
  (8 GB RAM, `-profile local` with the machine-local resource config): 24 of
  24 tasks, 3 min 12 s, including a Salmon index of the full mussel
  transcriptome. Harmonizing the output resolved all 425 records `exact`
  against `mgall_gene_id_crosswalk.tsv` under taxid 29158. All 60 simulated
  changes were called at padj < 0.05 in the right direction, 3 of 340
  unchanged genes were called, and estimated and true log2 fold changes
  correlated at 0.989. The simulated reads were not committed.
- **The manifest now checksums its reference.** The run showed the RNA-seq
  manifest named the assembly but recorded no checksum for the transcriptome,
  GTF or tx2gene it read. It now stages and hashes them, compares each with
  the bundle's `reference_provenance.json`, and warns on a mismatch, on a
  bundle for a different assembly, or on reference files with no bundle
  provenance.

Not done: bundles for the other five species are fetched on demand when a
study needs them (each `--plan` already resolves). The methylation manifest
does not yet checksum its genome FASTA, and the Bismark index is still built
by the workflow rather than kept in the bundle.

- Write a `data/reference/README.md` recipe, or an `aree fetch-reference
  --assembly GCF_...` command, that downloads the RefSeq transcript FASTA and
  GTF from the NCBI genomes FTP, verifies `md5checksums.txt`, builds
  `tx2gene.tsv`, and records a provenance JSON, generalizing the shell steps
  in [first_raw_reanalysis.md](first_raw_reanalysis.md). Salmon indexes are
  built by the workflow on first run.
- For the methylation workflow, extend the bundle with the genomic FASTA and
  a Bismark index, which the oyster bundle does not yet carry because no
  methylation run has completed.
- Storage: roughly 2 to 5 GB per species bundle before Salmon indexes; FASTQ
  volumes follow the chosen studies (the Manila clam deposits alone are
  several TB in aggregate, so stage per study as now).

Acceptance: a `config/<STUDY>.config` for one mussel RNA-seq study runs end
to end with `-profile local,workstation`, producing a manifest that names the
mussel assembly and annotation release.

### WP4. Second-species demo and regression fixtures

**Status: done, 2026-10-12.** `MYTGAL_HEAT01` is a simulated *M.
galloprovincialis* larval heat-stress study. It reports the oyster demo's gene
identifiers under the same phenotype (`larval_viability`) as the two oyster
larval studies, which pool *hsp70* and disagree on *sod1*, so without the
species partition it would have joined both pools. What was done and found:

- **The partition holds.** `tests/test_demo_species.py` runs the demo with and
  without the mussel study and requires every oyster meta-analysis row and
  every oyster candidate, down to adjusted p-values and scores, to be
  identical. They are. The mussel study forms its own pools (k = 1), and its
  candidates never reach the high-priority tier.
- **Seven existing tests assumed one species.** Each took the first row or card
  for a gene and now got the mussel one. No oyster result changed; the tests
  now select the oyster partition explicitly (`oyster()` and `oyster_card()` in
  `tests/conftest.py`).
- **Two outputs could not tell species apart.** Evidence cards printed only a
  taxid, and the top-candidates summary listed the oyster and mussel *hsp70*
  as two identical-looking rows. Cards now name the species and state that they
  pool one species only; the summary has a species column.
- **The demo list is no longer hard-coded.** It was spelled out in the
  Makefile, both GitHub workflows, the README and the test fixtures. It is now
  every `simulated: true` study (`aree study-ids --simulated`), so the next demo
  study needs no list edits.
- The demo study declares `species_identity` (collection site) and carries
  `species_assignment_unverified`, exercising the WP0 species-complex rule.
- The synthetic crosswalk header says both demo species share its identifiers
  on purpose, and drops its stale `AREE_CROSSWALK` instructions.

- Add one simulated *M. galloprovincialis* study (`MYTGAL_HEAT01`, clearly
  `simulated: true`) sharing feature identifiers with the oyster demo so the
  species partition is exercised by `make demo`, the dashboard, and the
  evidence-card tests. This is the cheapest way to keep multi-species
  behaviour from regressing before real data lands.
- Add a demo mussel entry to the synthetic crosswalk notes so nobody mistakes
  demo LOC numbers for real mussel genes, as the oyster demo already warns.

### WP5. First real studies (curation backlog)

**Status: backlog done, 2026-10-12; Manila clam heat pool done, 2026-10-17
(k = 2: TUCCI2025_HEAT_PRIMING, PERUZZA2023_HEATWAVE, both 3' Tag-seq,
fitted on clams with tank_level_replication). No mussel study yet.** The screen is in
[candidate_studies_mussels_clams.md](candidate_studies_mussels_clams.md): 186
ENA RNA-seq projects listed, about 40 inspected from deposited metadata, and
publications matched in PubMed. Three within-species pools are within reach:
Manila clam heat (k = 2, about 22 GB in total, curate first), hard clam
salinity (k = 3), and *M. edulis* heat (k = 2, up to 4). Findings:

- **Mixed assays inflate group sizes.** Both Manila clam heat projects deposit
  16S amplicon runs from the same animals under the same attributes, which
  doubled every group in `fetch-samplesheet`. It now takes
  `--library-strategy` and, whenever a project mixes strategies, prints group
  sizes per strategy instead of filtering silently (the label can be wrong,
  as it was for `LUTIER2022_OA_TIPPING`).
- **Three hard clam acidification projects look like one experiment** (one
  Stony Brook rearing program sampled at three ages), so they should be one
  study, not k = 3.
- **Pooling within *Mytilus* is thin.** *M. galloprovincialis* has no
  same-phenotype pair with n ≥ 3; its strongest study (larval acidification
  with measured abnormality) waits for cross-species pooling in WP6.
- The one candidate with a measured survival difference between sequenced
  groups is Manila clam heat priming (Tucci et al. 2025), primed vs naive.

Open `docs/candidate_studies_mussels_clams.md` using the same screen as the
oyster backlog: raw reads public, n ≥ 3 biological replicates confirmed from
run metadata, a measured resilience outcome preferred over exposure-only,
RNA-seq first, and stressor spread so pooled groups reach k ≥ 2. Verify every
accession against ENA before listing it. Candidate dataset classes, by
species and context, are in section 6. Each study is then claimed with
`aree start-study` on its own branch exactly as oyster studies are.

Acceptance: at least two mussel and two clam RNA-seq studies harmonized with
real standard errors, giving at least one within-species pooled group with
k ≥ 2 per taxon.

### WP6. Cross-species orthology and pooling (the new science)

**Status: done. Orthogroups 2026-10-13; pooling and cross-species candidates 2026-10-14; dashboard gate 2026-10-17.** What exists:

- `aree fetch-reference --proteins` writes one protein per gene (the longest)
  for each species, keyed by NCBI GeneID; all seven reference assemblies have
  one, 100% covered by their crosswalks (228,833 proteins).
- `aree build-orthogroups` ran OrthoFinder 3.1.5 on them (6 threads, 55
  minutes on Pori) in its own conda environment. The run,
  `orthofinder3.1.5_2026-10-13`, is committed and active: 38,507 orthogroups,
  of which 5,552 are single-copy in all seven species, 7,031 single-copy where
  present, 11,877 multi-copy and 14,047 species-specific. Shared with at least
  one other species: 95-98% of each mussel's genes, 85-88% of each clam's, 77%
  of the oyster's.
- Harmonization attaches the active run's orthogroup to every resolved real
  GeneID and records the run in its manifest. On HESSER2024_VCOR, 274 of 306
  resolved genes got one (the rest are non-coding), and 190 of those 274 sit
  in multi-copy orthogroups, only 18 in single-copy ones. The worry below
  about lineage-specific expansions is therefore the normal case, not the
  exception: most genes stress studies report are in expanded families.
- Checks: in 63% of single-copy-in-all orthogroups with named genes, every
  species' gene name agrees exactly (the rest are naming variants such as
  "forkhead box C1-a" vs "forkhead box protein C1"); "heat shock 70 kDa
  protein" genes fall into 285 orthogroups, the largest holding 17-23 copies
  in each mussel against 3 in the oyster.

**Orthogroup-level meta-analysis: done, 2026-10-14.** `aree meta-analyze
--level orthogroup` pools each species' gene-level estimate across species
through orthogroups (two-stage DerSimonian-Laird), reports between-species
Q, I-squared and tau-squared, per-species estimates and study counts, and a
`meets_replication_gate` flag (two or more species with k >= 2 each).
Multi-copy orthogroups are set aside by default, as are orthogroups missing
from the active run. On the demo, 8 orthogroups pool oyster with mussel
(larval viability), none past the gate; on all committed real evidence
(465,531 oyster records), nothing pooled across species yet, correctly, until
the first clam study merged (see below).
Building it found that the demo's illustrative orthogroup IDs collide with
real OrthoFinder IDs, so simulated evidence never consults the run (tested).

**Cross-species candidates and cards: done, 2026-10-14.** `aree
build-cross-species-cards` tiers each orthogroup row behind five named gates
(replication, direction, significance, phenotype relevance, quality) into
`cross_species_convergence`, `cross_species_emerging`,
`cross_species_discordant` or `cross_species_unsupported`, with a transparent
weighted score and a between-species I-squared penalty. Each card shows every
species' estimate side by side, the pooled estimate, the orthogroup run and
the gates failed. Single-species cards now name their orthogroup. On the
demo: 5 emerging, 3 unsupported, none convergent.

First real cross-species rows, 2026-10-16, after TUCCI2025_HEAT_PRIMING
merged: oyster and Manila clam thermal-tolerance evidence pools through 5,681
single-copy orthogroups, giving 64 emerging, 6 discordant and 5,611
unsupported rows, none convergent, with 603 cards. Real data changed two
rules. "Discordant" now requires species that are each significant at gene
level and disagree in sign: 2,880 of 2,886 sign disagreements involved a
non-significant species, which is noise, not a reversal. Quality is now
judged per species, because the union of flags across species always
exceeded the limit.

With PERUZZA2023_HEATWAVE merged (2026-10-17) the clam reached k = 2 and the
replication gate is met on real data for the first time. On all of `main`'s
real evidence, oyster and clam pool through 5,975 single-copy orthogroups
(5,152 multi-copy set aside): 133 emerging and 5,842 unsupported rows, none
discordant or convergent, with 455 cards. The gate is met in 5,633 of them;
131 emerging rows fail only phenotype relevance and quality. Both clam
comparisons are exposure-only, so no oyster-clam row can converge on
thermal tolerance until a clam study measures a tolerance outcome in its
pooled comparison. The 6 discordant rows seen with one clam study did not
survive the second.

**Dashboard gate: done, 2026-10-17.** The dashboard's pools are now per
species, and `cross_species_gate` opens a phenotype only when at least two
species each have a within-species gene pool with k >= 2. For an open
phenotype the home page's "Across species" section shows the gate, the
cross-species tiers and the top rows with their cards; rows for a closed
phenotype are counted as withheld, never ranked or linked. Whether or not the
gate is open, each top single-species candidate carries an ortholog lookup:
the other species' genes in its orthogroup, the phenotypes they were tested
under, and how many of their records were significant. On `main` (2026-10-17)
thermal tolerance is open (*M. gigas* k = 3, *R. philippinarum* k = 2) and
the demo's larval viability is closed (the demo mussel has k = 1, 8 rows
withheld). `make demo`, CI and the Pages build now run
`aree build-cross-species-cards`.

The
OrthoFinder run used `-M dendroblast` (OrthoFinder 3 has no orthogroups-only
mode) and only its `Orthogroups.tsv` is used; hierarchical orthogroups are a
possible refinement.

- **Orthogroup inference.** Run OrthoFinder (or use OrthoDB) on the RefSeq
  proteomes of every species with a committed crosswalk, starting from the
  oyster. Add `aree build-orthogroups` that writes
  `data/reference/orthogroups/<run>/` with the Orthogroups table, software
  version, proteome checksums, and a provenance JSON, then populates
  `orthogroup_id` in each species crosswalk from that run. The crosswalk's
  empty-orthogroup invariant becomes "empty unless an orthogroup run is
  recorded in provenance".
- **Mapping confidence gets its meaning back.** `one_to_one_ortholog`,
  `one_to_many_ortholog` and `many_to_one_ortholog` are currently used for
  within-species identifier ambiguity. Record the orthogroup's per-species
  copy number alongside, so a gene in a single-copy orthogroup across all
  pooled species is distinguishable from one in a 40-member expansion.
  Bivalve genomes are full of lineage-specific expansions (HSP70, C1q, IAPs,
  GIMAPs), and those are exactly the families stress studies report, so most
  "interesting" genes will be one-to-many. Pooling must default to
  single-copy orthogroups and must label the rest.
- **Meta-analysis at orthogroup level.** Add `aree meta-analyze --level
  orthogroup`. Group keys become `orthogroup_id + phenotype + feature_type +
  simulated`, with `species_taxid` retained as a subgroup so the output
  reports per-species pooled estimates, the cross-species estimate, and
  between-species heterogeneity separately. Effect sizes stay log2 fold
  changes; do not pool across assay types any more than today.
- **Prioritization.** Add a fourth candidate category, *cross-species
  convergence*: same orthogroup, same phenotype, consistent direction in at
  least two species, each species itself supported by at least one study
  with real standard errors. The score components already listed in
  [interpreting_candidate_scores.md](interpreting_candidate_scores.md) extend
  with number of species and orthogroup copy-number class. Evidence cards
  show per-species forest-plot rows and name the orthogroup run.
- **Gate.** No cross-species pooling is shown on the dashboard until at least
  two species each have a within-species group with k ≥ 2 on the same
  phenotype. Before that the orthogroup view is a lookup ("where is this
  oyster candidate's ortholog in mussels, and was it tested?"), not a
  meta-analysis.

### WP7. Identifier crosswalks for non-RefSeq genomes (second wave)

For *M. coruscus*, *S. constricta* and *P. generosa* there are no NCBI Gene
records. Add a second builder path, `aree build-crosswalk --from-gff
<annotation.gff3> --proteins <proteins.faa>`, that takes gene IDs from the
authors' annotation, attaches UniProt accessions through a recorded
DIAMOND/BLASTP run against UniProtKB with explicit thresholds, and marks every
such row `inferred` rather than `exact`. Orthogroups from WP6 then give these
species their only route to cross-species comparison. Until this lands, their
studies harmonize with `unresolved` identifiers and are visible as registered
but unmapped.

### WP8. Interface and documentation

- Dashboard and study browser: species filter, per-species crosswalk and
  study counts, a species-identity column for *Mytilus* studies.
- Docs: turn [adding_a_species.md](adding_a_species.md) into a worked example
  using the first mussel; add a page on species identity and
  reference-species mismatch; update [roadmap.md](roadmap.md),
  [implementation_status.md](implementation_status.md) and
  [identifier_mapping.md](identifier_mapping.md) (orthogroup section).

## 6. Dataset classes to curate first

Described by type and resilience context, not by accession; every accession
must be verified against ENA and the paper before it goes into the backlog
page, and the ENA counts in section 3.1 say roughly how much there is to
screen.

**Mussels**

1. *M. galloprovincialis* versus *M. trossulus* heat-stress gill RNA-seq from
   their Pacific contact zone (congeneric thermal-tolerance contrast). The
   best-known datasets are microarray, so look for later RNA-seq replications.
   This is the first cross-species pair WP6 would pool.
2. *M. edulis* / *M. galloprovincialis* mass-mortality or summer-mortality
   cohorts with survivors versus dead, or resistant versus susceptible
   families (Atlantic France, Mediterranean, Adriatic). A measured survival
   phenotype.
3. *M. galloprovincialis* warming × acidification multi-stressor transcriptomics
   (Mediterranean), paired with the oyster OA studies already registered so
   `acidification_tolerance` reaches k ≥ 2 within and then across species.
4. *M. californianus* intertidal heat and emersion RNA-seq, and its two
   bisulfite deposits (94 runs): the first mussel methylation study and the
   first use of `emersion_desiccation`.
5. Mussel hemocyte RNA-seq under *Vibrio* challenge, resistant versus
   susceptible lines where available (`immune_responsiveness`,
   `disease_resistance`).
6. *M. coruscus* hypoxia and heat selective-breeding lines (second wave;
   register now, harmonize after WP7).

**Clams**

7. *R. philippinarum* brown ring disease (*Vibrio tapetis*) resistant versus
   susceptible hemocyte or gill RNA-seq. Measured resistance phenotype.
8. *R. philippinarum* *Perkinsus olseni* infection load transcriptomics
   (`pathogen_load`).
9. *R. philippinarum* heat, hypoxia, air exposure and low-salinity RNA-seq,
   including orange/white or other breeding strains where strain survival was
   scored.
10. *M. mercenaria* QPX-resistant versus susceptible lines, and heat ×
    hypoxia × acidification multi-stressor RNA-seq (`multi_stressor`).
11. *Mya arenaria* heat or hypoxia transcriptomics; its transmissible
    neoplasia RNA-seq is disease-associated only.
12. *P. generosa* larval and juvenile OA transcriptomics and methylation from
    this lab (second wave for identifiers; check how the bisulfite runs are
    labelled in SRA, since none were found under Bisulfite-Seq).

Priority inside each list: the studies that give a within-species group k ≥ 2
on `thermal_tolerance`, `acidification_tolerance` or `disease_resistance`,
because those are the phenotypes where oyster evidence already exists and
WP6 could first be tested.

## 7. Sequencing and rough effort

| Order | Package | Depends on | Effort (person-days, rough) |
|---|---|---|---|
| 1 | WP0 vocabularies and reference data | none | 1 to 2 |
| 2 | WP1 per-study crosswalk selection | WP0 | 2 to 3 |
| 3 | WP2 six crosswalks | WP1 | 1, mostly download time |
| 4 | WP4 simulated mussel demo | WP1 | 1 |
| 5 | WP3 reference bundles + first config (done) | WP0 | 1 to 2 plus a workflow run |
| 6 | WP5 backlog page, then first four studies | WP1 to WP3 | 1 for the page; each study as today (a day or more of compute per raw reanalysis) |
| 7 | WP6 orthology and cross-species pooling | WP2, WP5 | 5 to 10, plus an OrthoFinder run (hours to a day for ten proteomes) |
| 8 | WP7 non-RefSeq crosswalks | WP6 | 3 to 5 |
| 9 | WP8 interface and docs | ongoing | 2 to 3 |

Milestone A (WP0 to WP4): mixed-species pool runs in CI with simulated data.
Milestone B (WP5): first real mussel and clam evidence on the dashboard,
species-partitioned. Milestone C (WP6): first cross-species convergence
candidate with an evidence card that shows both species' contributions and
the orthogroup copy-number class.

## 8. Assumptions and open decisions

- **Assumption:** NCBI GeneID remains the anchor identifier for every species
  on the RefSeq route, so the identifier hierarchy and retired-ID handling
  carry over unchanged. Ensembl Metazoa coverage for these species is thin
  and is not relied on.
- **Assumption:** The oyster RNA-seq and methylation workflows need no
  species-specific logic beyond reference paths; mussel and clam libraries
  are ordinary poly-A or rRNA-depleted Illumina data. Tag-seq and
  microarray deposits are excluded as for oysters.
- **Decision needed:** whether geoduck should jump the queue ahead of the
  RefSeq species despite needing WP7, given lab priorities. This plan says no,
  but registers geoduck studies early so they are counted.
- **Decision needed:** whether the three *M. edulis*-complex species should
  ever be pooled below orthogroup level (they share most GeneID-level
  annotation structure but not GeneIDs). This plan says no: pooling across
  them goes through WP6 like any other pair, with the complex shown as a
  display grouping only.
- **Decision needed:** compute and storage for the first-wave FASTQ staging,
  which is larger than the oyster backlog because of the Manila clam volume.
- **Risk:** species mislabelling in *Mytilus* deposits is likely common
  enough that some first-wave studies will carry
  `species_assignment_unverified` permanently. That is acceptable as long as
  it is visible; it is not acceptable to drop the flag to make a pool larger.
- **Risk:** lineage-specific gene-family expansions mean the cross-species
  pool of single-copy orthogroups may exclude most of the genes stress
  studies headline. That is a true statement about the biology, and the
  evidence cards should say so rather than pool paralogs.
