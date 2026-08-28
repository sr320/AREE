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

AREE ships **two** crosswalks and never merges them.

| File | Contents | When it is used |
|---|---|---|
| `data/mappings/gene_id_crosswalk.tsv` | 21 synthetic genes backing `data/demo/*` | Default. Keeps `make demo` and the test suite self-contained. |
| `data/reference/crosswalk/mgigas_gene_id_crosswalk.tsv` | 33,356 real *Magallana gigas* genes from NCBI Gene + UniProtKB | Curating real studies. Selected explicitly. |

Selection precedence is `set_crosswalk_path()` > `$AREE_CROSSWALK` > the demo
crosswalk. To harmonize real data:

```bash
aree build-crosswalk                       # writes the real crosswalk + provenance sidecar
export AREE_CROSSWALK=data/reference/crosswalk/mgigas_gene_id_crosswalk.tsv
aree harmonize --study <STUDY_ID> --input <results.tsv>
```

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
  Ensembl xrefs, description, gene type.
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
3. **`orthogroup_id` is empty by design.** Real orthogroups require an ortholog
   inference run (OrthoFinder/OrthoDB). Fabricating them would manufacture
   cross-species confidence that does not exist, so the column ships empty and
   `resolve_identifier` returns `None` for it.

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

- **Stable ID columns** (`ncbi_gene_id`, `ensembl_gene_id`,
  `uniprot_accession`, `locus_id`): looked up directly against the crosswalk.
  A single match is `exact`; multiple matching rows collapse to
  `many_to_one_ortholog`; no match is `unresolved`. Two of these get extra
  handling for real data: `ensembl_gene_id` matches within a `;`-delimited list,
  because NCBI records more than one Ensembl xref for some genes, and
  `uniprot_accession` is resolved through a prebuilt index over
  `uniprot_accessions_all` so that any accession for a gene — not just the
  representative one stored in `uniprot_accession` — resolves `exact`.
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

## Related documentation

- [interpreting_candidate_scores.md](interpreting_candidate_scores.md)
- [adding_a_species.md](adding_a_species.md) — crosswalk files are per-species
- [design.md](design.md#5-identifier-harmonization)
