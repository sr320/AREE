# Identifier mapping

Different assay types and different studies report molecular features under
different identifier systems — NCBI Gene IDs, Ensembl IDs, UniProt
accessions, legacy locus IDs, or bare gene symbols. AREE never overwrites the
original identifier; it adds a standardized identifier alongside it, with an
explicit confidence label for how much to trust the translation.

## The identifier hierarchy

From highest to lowest precedence, implemented in
`src/harmonize/identifiers.py`:

```
NCBI Gene ID > Ensembl gene ID > UniProt accession > locus ID (reference oyster annotation) > gene symbol > orthogroup
```

Every evidence record keeps both:

- `feature_id_original` — exactly as reported in the source result file, never
  altered.
- `feature_id_standardized` — the harmonized identifier resolved via this
  hierarchy, or `null` if `mapping_confidence` is `unresolved`.

## Which crosswalk is in force

AREE keeps a synthetic demo crosswalk and one real crosswalk per species, and
never merges any of them.

| File | Contents | When it is used |
|---|---|---|
| `data/mappings/gene_id_crosswalk.tsv` | 21 synthetic genes backing `data/demo/*` | Every study with `simulated: true`. Keeps `make demo` and the test suite self-contained. |
| `data/reference/crosswalk/mgigas_gene_id_crosswalk.tsv` | 33,356 real *Magallana gigas* genes from NCBI Gene + UniProtKB | Real studies on either oyster assembly. |
| `data/reference/crosswalk/{mgall,medul,mtros,mcali,mmerc,rphil}_gene_id_crosswalk.tsv` | Real genes for six mussel and clam species, same sources and builder (see below) | Real studies on that species' assemblies. |

**The study decides, through its genome assembly.** Each entry in
`data/reference/genome_assemblies.yaml` names its `identifier_crosswalk`.
Harmonizing a real study looks up the study's `genome_assembly` there and
resolves against that file, so an oyster study and a mussel study harmonized in
the same run each use their own species' identifiers. No environment variable is
needed:

```bash
aree harmonize --study <STUDY_ID> --input <results.tsv>
```

The assembly decides rather than the study's species because identifiers come
from the annotation the reads were mapped to. A study that mapped reads to a
related species' genome (quality flag `reference_species_mismatch`) therefore
reports the reference species' gene identifiers, while its evidence is still
partitioned under the study animal's `species_taxid`.

A real study whose assembly has no registered crosswalk is an error, not a
fallback to some other crosswalk. The error names the `aree build-crosswalk
--taxid <n>` command to run, after which the new file is registered by setting
`identifier_crosswalk` for that species' assemblies.

`$AREE_CROSSWALK` still exists as an override for local experiments, such as
testing a rebuilt crosswalk before committing it. It is honoured only if its
provenance sidecar shows it was built for the study's assembly taxon; a
crosswalk for another species, or one with no provenance, stops harmonization.
Earlier versions of AREE let the most recently selected crosswalk carry over to
the next study in the same process. That is no longer possible: every study
selects afresh.

**The two files must never be concatenated.** The demo crosswalk's LOC numbers
are not unused placeholders — most of them collide with real NCBI GeneIDs that
denote completely different genes. `LOC105333935` is labelled `hsp70` in the
demo file but is really *ubiquitin carboxyl-terminal hydrolase 32*;
`LOC105341000` is labelled `actin` but is really *ADAMTS adt-1*. A union would
silently attach demo biology to real accessions, and nothing downstream would
flag it.

## The real crosswalk

Built by `src/mappings/build_crosswalk.py` (`aree build-crosswalk`) from two
public sources, with a `.provenance.json` sidecar recording source URLs, SHA-256
checksums of both inputs and the output, retrieval date, and coverage:

- **NCBI Gene `gene_info`**, filtered to taxid 29159 — GeneID, symbol, synonyms,
  Ensembl xrefs, description, gene type, and `LocusTag` (written as
  `legacy_locus_tags`; see [Legacy gene models](#legacy-gene-models-oyster_v9-cgi_-ids)).
- **UniProtKB REST**, same organism — accession, reviewed status, GeneID xref.
- **NCBI Gene `gene_history`**, same organism — GeneIDs NCBI has discontinued and
  what replaced them (written to a separate sidecar, see below).

Coverage as built (2026-08-27, NCBI annotation release `GCF_963853765.1-RS_2024_06`):

| Column | Genes covered | % |
|---|---|---|
| `ncbi_gene_id`, `locus_id` | 33,356 | 100% |
| `ensembl_gene_id` | 12,717 | 38.1% |
| `uniprot_accession` | 2,790 | 8.4% |
| named `gene_symbol` (non-LOC) | 650 | 2.0% |
| `legacy_locus_tags` (oyster_v9 `CGI_` models; added in builder v1.1.0, 2026-10-17) | 14,699 | 44.1% |
| `orthogroup_id` | 0 | 0% |

Plus a sidecar, `mgigas_retired_gene_ids.tsv`, covering 9,057 discontinued
GeneIDs that NCBI has remapped to a current gene.

Three of these numbers deserve attention before anyone reads a proteomics
result:

1. **UniProt linkage is sparse (8.4%).** Of 58,641 UniProtKB entries for this
   organism, only 3,170 carry any NCBI GeneID cross-reference; 55,501 carry
   neither a GeneID nor a RefSeq xref. This is a property of the source data,
   not of AREE. Proteomics evidence will therefore have a materially higher
   `unresolved` rate than transcriptomics evidence, and that difference must not
   be read as a biological signal.
2. **Only 2% of genes have a real name.** The rest are `LOC`-form. Symbol-based
   matching is close to useless for this species; identifier-based matching is
   the only reliable path.
3. **`orthogroup_id` is empty in every real crosswalk.** Crosswalks hold only
   what NCBI and UniProt provide. Orthogroups come from a separate inference
   run (see *Orthogroups* below); harmonization attaches them to evidence from
   the active run, never from the crosswalk.

### Orthogroups

Orthogroups link genes across species. They are inferred, not looked up, so
they live apart from the crosswalks: `aree build-orthogroups` runs OrthoFinder
on one protein per gene for every species with a crosswalk (the
`primary_proteins.faa` that `aree fetch-reference --proteins` writes; see
[reference_bundles.md](reference_bundles.md)). Because those proteomes come
from the same RefSeq annotations the crosswalks are built on, every orthogroup
member is a GeneID the harmonizer already resolves.

A run is committed under `data/reference/orthogroups/<run_id>/`:

| File | Contents |
|---|---|
| `gene_orthogroups.tsv` | `ncbi_gene_id`, `ncbi_taxid`, `orthogroup_id` for every gene |
| `orthogroups.tsv` | per orthogroup: `copy_class`, species and gene counts, copies per species |
| `provenance.json` | OrthoFinder version and command, each proteome's assembly and SHA-256, per-species assignment rates |

`data/reference/orthogroups/active.yaml` names the run harmonization uses.
`resolve_identifier` attaches that run's orthogroup to every resolved real
GeneID, and each harmonization manifest records the run ID and the checksum of
its gene table. The synthetic demo crosswalk keeps its own illustrative
orthogroups and is never matched against a real run.

`copy_class` says how comparable an orthogroup is across species:

| Class | Meaning |
|---|---|
| `single_copy_all` | exactly one gene in every species of the run |
| `single_copy` | at most one gene per species, in two or more species |
| `multi_copy` | two or more genes in at least one species |
| `species_specific` | genes from one species only |

Bivalve genomes carry large lineage-specific expansions (HSP70, C1q, IAPs), and
those are the families stress studies report most. A `multi_copy` orthogroup
links a gene to all of its paralogs in another species, not to one
counterpart, so cross-species comparisons must say which class they rest on.

The active run, `orthofinder3.1.5_2026-10-13` (OrthoFinder 3.1.5, 7 species,
228,833 proteins):

| Copy class | Orthogroups |
|---|---|
| `single_copy_all` | 5,552 |
| `single_copy` | 7,031 |
| `multi_copy` | 11,877 |
| `species_specific` | 14,047 |

On a real oyster study (HESSER2024_VCOR), 274 of 306 resolved genes received
an orthogroup (the rest are non-coding genes, which have no protein), and 190
of those 274 are in `multi_copy` orthogroups. Expect most cross-species links
for stress-response genes to be one-to-many.

### Mussel and clam crosswalks

Built 2026-10-11 with the same builder, each against its species' NCBI
reference annotation (recorded in the provenance sidecar and checked by a
test). NCBI `gene_info` and `gene_history` were downloaded once and filtered
per taxid; every species' rows were confirmed contiguous in both files, which
the builder's streaming shortcut assumes.

| Crosswalk | Species | Genes | Ensembl xref | Named symbol | UniProt-linked genes | Retired GeneIDs remapped / dead |
|---|---|---|---|---|---|---|
| `mgall` | *Mytilus galloprovincialis* | 43,869 | 59.9% | 10.9% | 42 | 0 / 0 |
| `medul` | *Mytilus edulis* | 49,210 | 51.0% | 8.1% | 31 | 0 / 12 |
| `mtros` | *Mytilus trossulus* | 47,781 | 58.7% | 12.1% | 94 | 0 / 0 |
| `mcali` | *Mytilus californianus* | 41,309 | 55.0% | 9.4% | 7 | 0 / 0 |
| `mmerc` | *Mercenaria mercenaria* | 45,378 | 56.2% | 2.4% | 2 | 3,656 / 9,078 |
| `rphil` | *Ruditapes philippinarum* | 48,038 | 57.1% | 2.5% | 18 | 0 / 12 |

What this means in practice:

- **UniProt is effectively unusable for these species.** UniProtKB holds 81,997
  entries for *M. galloprovincialis* and 74,071 for *M. edulis*, but only 43 and
  32 of them cross-reference an NCBI GeneID. They are translated from GenBank
  submissions, not from the RefSeq annotation the crosswalk is built on. A
  mussel or clam proteomics table reported in UniProt accessions will resolve
  almost nothing until a sequence-based mapping exists. Report RefSeq protein
  or GeneID identifiers where the source allows.
- **Ensembl xrefs are unusually good** (51 to 60%, against 38% for the oyster),
  so a study reporting Ensembl gene IDs resolves well.
- **Hard clam has two annotated assemblies.** Both use `mmerc`. Of the 43,960
  GeneIDs on the older assembly (ASM1480567v1.1), 3,656 were discontinued in
  February 2023 and remap through `mmerc_retired_gene_ids.tsv`, and 9,078 have
  no replacement. See that assembly's notes in
  `data/reference/genome_assemblies.yaml`.
- **No retirement history yet** for most mussels: their RefSeq annotations are
  first releases, so there are no older GeneIDs to remap.

### Genomes without NCBI Gene records (GFF route, WP7)

A species NCBI has not annotated has no GeneIDs, so the NCBI route has nothing
to download. `aree build-crosswalk --from-gff` builds its crosswalk from the
authors' own annotation instead (`src/mappings/gff_crosswalk.py`):

- **Genes come from the authors' GFF3**, one row each, anchored on `gene_id`,
  the GFF gene ID without a `gene-` prefix. For an NCBI GenBank annotation that
  is the locus tag (`MCOR_1`). `ncbi_gene_id` stays empty. Every downstream
  reader uses the crosswalk's anchor, never a route's own column, and the
  provenance records `"route": "gff_homology"` and `"anchor_column": "gene_id"`.
- **Transcript and protein IDs** from the GFF go into `transcript_ids` and
  `protein_ids` (GenBank protein IDs with and without their version). A study
  reporting any of them reaches the gene.
- **UniProt accessions come from one recorded DIAMOND blastp search** of the
  longest protein per gene against a named UniProtKB release (Swiss-Prot by
  default). The thresholds are explicit and recorded: e-value ≤ 1e-10,
  identity ≥ 30%, query and subject coverage ≥ 50%, up to 5 hits per gene. The
  builder re-applies them when it reads the hits table. The best hit is
  `uniprot_accession`, all passing hits are in `uniprot_accessions_all`, and
  `uniprot_evidence` is `diamond_blastp`. The source files, the database and
  the hits table are all checksummed in the provenance.
- **The authors' own homology names** (GenBank `standard_name`) become
  `gene_symbol`, and their own UniProt assignment, where the GFF records one,
  is kept in `author_uniprot_accession`. The provenance reports how often it
  is among our hits, a free check on the search.

Resolution follows what a lookup went through, not the route that built the
file:

| Study reports | Resolves |
|---|---|
| the authors' gene ID, locus tag, transcript ID or protein ID | `exact`: it names that gene in that annotation |
| a UniProt accession | `inferred`: reached through our homology search; `many_to_one_ortholog` when several genes share it |
| a gene symbol | `inferred` or `many_to_one_ortholog`, as for every route |

A GFF-route gene gets no orthogroup from the active orthogroup run, which is
keyed by NCBI GeneID. An author ID could collide with a GeneID, so the lookup
is skipped. Cross-species comparison needs the species added to an orthogroup
run, which is not built yet.

**First GFF-route crosswalk: `mcoru`**, *Mytilus coruscus* MCOR1.1
(GCA_011752425.2, GenBank annotation submitted by Ningbo University). Built
2026-10-17 against UniProtKB/Swiss-Prot 2026_03 with DIAMOND 2.2.8
(`--sensitive`, 103 s on 6 threads).

| Genes | With a protein | With a Swiss-Prot hit | Authors' accession among our hits |
|---|---|---|---|
| 58,540 | 58,249 | 9,877 (17.0%) | 3,085 of the 6,353 genes both assigned (48.6%) |

The hit rate is low because the thresholds require half of both proteins to
align, which keeps partial, single-domain matches out. Resolving a study's own
gene, transcript or protein IDs does not depend on it. Of the 3,268 genes
where the authors' accession is not among our hits, 1,883 of the 3,265 that
can be compared (57.7%) carry the same Swiss-Prot protein name on both
entries: the same protein, recorded from a different species. Counting those,
our search names the authors' protein for 78% of the genes both assigned. *M. coruscus* studies on the two later chromosome-level
assemblies use annotations NCBI does not hold, so their identifiers will not
resolve against `mcoru`.

### Design decisions in the real crosswalk

- **One row per NCBI GeneID.** GeneID anchors the table because it is the most
  stable identifier across annotation releases. Emitting one row per
  (gene, protein) pair would make a plain GeneID lookup return several rows,
  which the resolver would correctly but misleadingly downgrade to
  `many_to_one_ortholog`.
- **UniProt is one-to-many, so two columns carry it.** `uniprot_accession` holds
  a single representative (reviewed preferred); `uniprot_accessions_all` holds
  the full `;`-delimited set. The resolver indexes *every* accession, so a study
  reporting a non-representative accession still resolves `exact`. 302 genes
  have more than one accession; 14 accessions are attributed to more than one
  gene and are downgraded to `many_to_one_ortholog`.
- **`locus_id` is `LOC<GeneID>` for every gene**, including genes that now carry
  an official symbol. The LOC form unambiguously denotes that GeneID by NCBI
  convention, and studies analysed against an older annotation routinely report
  the LOC form for a gene that has since been named — `LOC105331241` is now
  *coadhesin*. This is the single largest source of recoverable matches across
  annotation versions (see [handling_genome_versions.md](handling_genome_versions.md)).
- **`gene_synonyms` is emitted but not yet consulted** by `resolve_identifier`.
  It is there for curators and for a future synonym-resolution step.

### Retired GeneIDs

*M. gigas* has been re-annotated several times, and a study analysed against an
older annotation reports GeneIDs that no longer exist in the current one. NCBI
retired 22,658 GeneIDs for this organism, and without handling them every one
would resolve as `unresolved` — silently discarding real evidence from exactly
the older studies AREE most wants to reanalyze.

`mgigas_retired_gene_ids.tsv` records what NCBI says happened to each:

| Outcome | Count | Resolution |
|---|---|---|
| Replaced by a current gene | 9,057 | resolves to the replacement, `inferred` |
| Discontinued with no replacement | 13,601 | stays `unresolved` |

`LOC105320012`, for instance, was discontinued and replaced by `105343733`; it
now resolves, where before it did not.

The confidence label is deliberately **`inferred`, not `exact`**, even though the
replacement comes from NCBI itself and is unambiguous. The remap crosses an
annotation version, and AREE's convention is that only a direct hit on the
current annotation earns `exact`. If you disagree with that call, it is one
constant in `src/harmonize/identifiers.py` — but note that changing it would
make version-remapped evidence indistinguishable from directly matched evidence
in every downstream meta-analysis.

Retired IDs whose replacement is itself absent from the crosswalk are dropped at
build time rather than chained, so no identifier resolves through a gene that
does not exist.

### Legacy gene models (oyster_v9 `CGI_` IDs)

Studies from roughly 2012–2016 report the gene models of the first *C. gigas*
assembly (oyster_v9, Zhang et al. 2012), named `CGI_10001851` and so on. Most
proteomics studies of that period searched spectra against that proteome, so
their protein tables are keyed by these IDs, sometimes as ProteinProphet groups
(`CGI_10000041,CGI_10008574`). They are not GeneIDs, and the retired-GeneID
sidecar does not cover them.

NCBI keeps them, though: a current gene whose model descends from an oyster_v9
model carries that model's ID in `gene_info`'s `LocusTag` field. Builder v1.1.0
writes that field to `legacy_locus_tags`. In the current annotation 14,699
genes carry one, and no tag is shared between genes. `resolve_identifier(raw,
"legacy_locus_tag")` looks the ID up there:

- a tag on one current gene resolves to that gene, `inferred` (it crosses an
  annotation version, as with retired GeneIDs);
- a protein group resolves only if every member names the same gene;
  otherwise it is `unresolved`, because a peptide-level group cannot be split;
- a model NCBI did not carry forward has no tag and stays `unresolved`.

`src/harmonize/proteomics.py` sends any `CGI_` ID (or group) down this route and
everything else down the UniProt one. For the per-oyster NSAF table of PRIDE
`PXD000835` (2014), 1,107 of 1,616 proteins (68.5%) resolve this way; the rest
stay `unresolved`, visibly, rather than being matched by name.

### Rebuilding

The build streams ~230 MB from NCBI and takes roughly 15–20 minutes on a typical
connection. The built crosswalk and its provenance sidecar are committed, so a
rebuild is only needed when the NCBI annotation is updated. Raw downloaded
sources land in `data/reference/crosswalk/_sources/` and are gitignored.

## `mapping_confidence` levels

Defined in `registry/controlled_vocabularies/mapping_confidence.yaml`, ranked
strongest to weakest:

| rank | id | meaning |
|---|---|---|
| 1 | `exact` | Direct match on a stable, unambiguous identifier |
| 2 | `one_to_one_ortholog` | Unambiguous one-to-one orthology relationship |
| 3 | `one_to_many_ortholog` | Source identifier maps to more than one candidate reference identifier |
| 4 | `many_to_one_ortholog` | Multiple source identifiers collapse onto a single reference identifier |
| 5 | `inferred` | Mapping derived indirectly (e.g. symbol match below ortholog-calling thresholds) |
| 6 | `unresolved` | No confident mapping found |

`unresolved` is a valid, expected outcome, not an error. Unresolved records
stay in the evidence table (nothing is dropped) but are excluded from
identifier-level pooling in meta-analysis (`src/meta_analysis/run.py` filters
`mapping_confidence != "unresolved"` before grouping), because there is no
stable identity to group on.

## How resolution actually works (`resolve_identifier`)

`resolve_identifier(raw_id, id_type)` in `src/harmonize/identifiers.py`
handles two families of input:

- **Legacy gene models** (`legacy_locus_tag`): oyster_v9 `CGI_` IDs, matched
  against `legacy_locus_tags` and labelled `inferred`; see
  [Legacy gene models](#legacy-gene-models-oyster_v9-cgi_-ids).
- **Stable ID columns** (`ncbi_gene_id`, `ensembl_gene_id`,
  `uniprot_accession`, `locus_id`): looked up directly against the crosswalk.
  A single match is `exact`; multiple matching rows collapse to
  `many_to_one_ortholog`; no match is `unresolved`. Two of these get extra
  handling for real data: `ensembl_gene_id` matches within a `;`-delimited list,
  because NCBI records more than one Ensembl xref for some genes, and
  `uniprot_accession` is resolved through a prebuilt index over
  `uniprot_accessions_all` so that any accession for a gene — not just the
  representative one stored in `uniprot_accession` — resolves `exact`.
- **GFF-route crosswalks** add the authors' `gene_id`, `transcript_ids` and
  `protein_ids` to the stable-ID search (all `exact`), and a UniProt match on
  a row whose `uniprot_evidence` is `diamond_blastp` is `inferred`; see
  [the GFF route](#genomes-without-ncbi-gene-records-gff-route-wp7).
- **Gene symbols**: looked up case-insensitively against the crosswalk's
  `gene_symbol` column first. A single match is labeled `inferred` (not
  `exact`) — symbols are not guaranteed unique or stable, so even a clean
  match gets a lower confidence than a stable-ID match. Multiple matches
  collapse to `many_to_one_ortholog`. If no crosswalk match is found at all,
  the symbol is checked against a curated exceptions file,
  `data/mappings/ambiguous_symbol_map.yaml`, before finally falling back to
  `unresolved`.

## Worked example: `hsp70` / `LOC105333935`

The crosswalk row for the canonical hsp70 locus (`data/mappings/gene_id_crosswalk.tsv`):

```
ncbi_gene_id   ensembl_gene_id     uniprot_accession  locus_id       gene_symbol  orthogroup_id
LOC105333935   ENSCGRG00000001     K1QN12             LOC105333935   hsp70        OG0000001
```

Three ways a source study might report this feature, and what AREE resolves:

1. **Source reports `ncbi_gene_id = LOC105333935` directly.** Stable-ID
   lookup finds exactly one row → `feature_id_standardized = LOC105333935`,
   `mapping_confidence = exact`.
2. **Source reports the plain symbol `hsp70`.** Symbol lookup finds exactly
   one row → `feature_id_standardized = LOC105333935`,
   `mapping_confidence = inferred` (a symbol match, not a stable-ID match,
   even though it happens to be unambiguous in this crosswalk).
3. **Source reports the symbol `hsp70-like2`** (as used by
   `GIGAS_LARV05`'s legacy gene-symbol-only annotation — see
   `registry/studies/GIGAS_LARV05.yaml`, which deliberately exercises this
   path). This does not match any `gene_symbol` value in the crosswalk
   directly, so resolution falls through to
   `data/mappings/ambiguous_symbol_map.yaml`:

   ```yaml
   - source_symbol: "hsp70-like2"
     standardized_id: "LOC105333935"
     mapping_confidence: "one_to_many_ortholog"
     note: "Symbol suggests a paralog of the canonical hsp70 locus; multiple candidate paralogs exist, best-guess target shown."
   ```

   Result: `feature_id_standardized = LOC105333935`,
   `mapping_confidence = one_to_many_ortholog`. The evidence record is kept
   and still contributes to the standardized `LOC105333935` group, but with a
   visibly lower confidence than case 1 — and `mapping_confidence_score` in
   candidate scoring (`src/prioritize/scoring.py`) weights it at 0.5 instead
   of 1.0.

The same ambiguous-map file also documents `cgi-actin` →
`many_to_one_ortholog` (multiple source actin paralogs collapsing onto one
reference locus) and `igf-like` → `inferred` (a fuzzy partial match to
`igf1`), which are the same mechanism applied to different studies' legacy
annotations.

## Why this matters for scoring and meta-analysis

- Meta-analysis groups evidence by `(feature_id_standardized, phenotype,
  feature_type)` — a `many_to_one_ortholog` or `one_to_many_ortholog` mapping
  means two source symbols may end up pooled as "the same" feature. This is
  the deliberate trade-off documented above, and it is visible in the
  candidate's `mapping_confidences` field, not hidden.
- `mapping_confidence_score` in candidate scoring takes the **worst**
  mapping confidence among all contributing records for a candidate
  (`src/prioritize/scoring.py`, `worst_mapping = min(...)`), so a candidate
  is only as trustworthy, identifier-wise, as its weakest link.

## Metabolites are linked, not mapped

Metabolite evidence keeps the reported metabolite name as its
`feature_id_standardized`; it is never forced into the gene hierarchy above.
Its connection to genes is a separate, weaker relationship (metabolite →
KEGG compound → reaction → KO → crg GeneID for enzymes; metabolite → compound
→ curated transporter KO in KEGG BRITE ko02000 → crg GeneID for transporters)
resolved at ranking time from
`data/mappings/metabolite_gene_links_*.tsv`, chosen per partition exactly as
the crosswalk is: the simulated demo table for simulated evidence, the table
registered under `metabolite_gene_links` in `genome_assemblies.yaml` for real
evidence. See [metabolite_gene_links.md](metabolite_gene_links.md).

## Related documentation

- [interpreting_candidate_scores.md](interpreting_candidate_scores.md)
- [adding_a_species.md](adding_a_species.md) — crosswalk files are per-species
- [design.md](design.md#5-identifier-harmonization)
