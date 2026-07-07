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

The crosswalk table backing this resolution is
`data/mappings/gene_id_crosswalk.tsv`. In the demo repository this crosswalk
is entirely **synthetic** — its header comment says so explicitly: identifiers
are internally self-consistent across `data/demo/*` but do not correspond to
real NCBI/Ensembl/UniProt records. A real deployment replaces this file with
an actual crosswalk built from the reference genome annotation (see
[handling_genome_versions.md](handling_genome_versions.md)).

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
  `many_to_one_ortholog`; no match is `unresolved`.
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
