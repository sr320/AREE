# Design note: linking metabolite evidence to genes

Status: **implemented (MVP).** Decisions in section 7 were confirmed on
2026-10-16.

## 1. Problem

The `multi_omics_convergence` tier (`src/prioritize/rank.py`,
`EvidenceIndex.supporting_layers` / `other_layers`) links molecular layers only
by an identical `feature_id_standardized`. RNA-seq, methylation and proteomics
evidence all standardize to an NCBI GeneID, so they can meet. Metabolite
evidence (`src/harmonize/metabolomics.py`) standardizes to the putative
metabolite *name* (`lactate`, `succinate`, ...), so it never shares an ID with
gene-level evidence. As a result, metabolomics can never take part in
convergence, whatever the data show.

Metabolites are not gene products, so a shared ID would be the wrong fix. The
honest link is biochemical: *this metabolite is a substrate or product of a
reaction that this gene's product catalyses*. That is a weaker claim than "the
same gene was measured twice". It must be labelled weaker, scored as weaker,
and explained on the card.

## 2. The link chain

```
metabolite name (as reported, preserved)
  └─ curated name → compound   (KEGG C-number, ChEBI ID)         [curation step]
       └─ compound → reaction   (KEGG REACTION, substrate/product side)
            └─ reaction → KO    (KEGG Orthology)
                 └─ KO → crg gene (KEGG GENES, organism code `crg`)
                      └─ crg gene = NCBI GeneID → AREE crosswalk / orthogroup
```

Facts checked on 2026-10-16 against the KEGG REST API:

- `crg` (T03920, *Magallana gigas*) has 26,085 proteins, and 10,575 of them
  have a KO. `link/ko/crg` returns 10,878 gene→KO pairs.
- `crg` gene IDs are NCBI GeneIDs. **10,874 of the 10,878 are current GeneIDs
  in AREE's RS_2024_06 crosswalk** (`data/reference/crosswalk/mgigas_gene_id_crosswalk.tsv`).
  None are retired, and 4 are unknown. The builder still routes every ID through
  the retired-ID sidecar and records the 4 that don't resolve. It does not drop
  them silently.
- Currency metabolites fan out enormously. For example, ATP (C00002) appears in
  810 KEGG reactions. Without a hub filter, ATP alone would "link" most of the
  metabolic genome to any ATP change.

## 3. Proposed files (`data/mappings/`)

| File | Contents |
|---|---|
| `metabolite_compound_map.tsv` | Curated `metabolite_name` (normalized, lowercase), `synonyms`, `kegg_compound`, `chebi_id`, `match_basis` (`exact_name` / `synonym` / `curator_judgement`), `curator`, `date`, `note`. Hand-curated and small. It holds only the names that appear in registered studies. |
| `currency_metabolites.yaml` | The explicit list of compounds that are never used as a link (ATP, ADP, AMP, NAD⁺/NADH, NADP⁺/NADPH, H₂O, H⁺, CO₂, O₂, Pi, PPi, CoA, NH₃), with the reason. These metabolites still appear as evidence. They just don't link. |
| `metabolite_gene_links_mgigas.tsv` | Built locally, **not committed** (decision A), one row per (compound, reaction, KO, gene): `kegg_compound`, `kegg_reaction`, `equation_side` (`left`/`right`), `ko`, `ncbi_gene_id`, `link_type`. The per-compound gene count and `link_confidence` are derived when the table is loaded, so they cannot drift from the rows. |
| `metabolite_gene_links_mgigas.provenance.json` | Same shape as the crosswalk provenance. It records the KEGG release (`info/kegg`), the access date, every REST URL fetched with its sha256, the builder version, the sha256 of the output, coverage counts, unresolved GeneIDs, and limitations. |
| `metabolite_gene_links_demo.tsv` (+ provenance) | **Simulated**, against the fabricated demo LOC IDs. It obeys the same rule as the demo crosswalk: it is never merged with the real table. |

A builder, `src/mappings/build_metabolite_links.py`, fetches **only** the
compounds listed in `metabolite_compound_map.tsv`. It makes one
`link/reaction/<cpd>` call per compound, then the reaction→KO and
reaction-equation lookups, plus a single `link/ko/crg` call. It never mirrors
KEGG. `genome_assemblies.yaml` gains a `metabolite_gene_links:` key next to
`identifier_crosswalk:`, so a real study resolves links through its assembly
and a simulated one through the demo table, exactly as identifiers do today.

**No change to the evidence schema.** Metabolite records keep
`feature_id_original` (the source feature ID) and `feature_id_standardized`
(the name). The link is resolved at ranking time from the table, so the link
table can be rebuilt without re-harmonizing any study.

## 4. Mapping confidence

The existing `mapping_confidence` vocabulary is reused for `link_confidence`.
A metabolite→gene link is **never** `exact` or `one_to_one_ortholog`:

| Situation | `link_confidence` |
|---|---|
| Compound reaches exactly one crg gene, through a KO transferred from KEGG's ortholog assignment | `inferred` |
| Compound reaches 2 or more crg genes (the usual case: paralogs, or several reactions) | `one_to_many_ortholog` |
| Name→compound match was `curator_judgement`, not an exact name or synonym | downgraded to `inferred` at most, and flagged |
| Compound reaches more than `MAX_GENES_PER_COMPOUND` genes (proposed: 25) | link **not credited**. Recorded as `hub_compound` |
| Compound in `currency_metabolites.yaml` | link **not credited** |

Only metabolites at MSI annotation level 1–2 have a standardized name today
(level 3–4 → `unresolved`), so only those can link. That gate stays as it is.

`link_type` is `enzyme` for MVP. Transporter links (for example, KEGG BRITE
ko02000 / TCDB) are listed as future work. KEGG REACTION rarely encodes
transport, and inventing a transporter link would overstate what is known.

## 5. Ranking changes (`src/prioritize/rank.py`)

A second, **separately named** kind of convergence. It never merges into the
first.

- `EvidenceIndex` loads the partition's link table and builds both directions:
  gene → linked metabolite names, and metabolite name → linked GeneIDs.
- New per-candidate fields:
  - `linked_layers`: layers with a significant record (adj. p ≤ 0.05) under the
    **same phenotype and partition** for a *linked* feature, excluding layers
    already in `supporting_layers`.
  - `n_linked_layers`, `linked_features` (for example,
    `lactate→LOC… via R00703/K00016`), `link_confidences`,
    `link_table_sha256`.
  - `convergence_basis`: `shared_identifier`, `metabolite_gene_link`, or `none`.
- New tier `metabolite_linked_convergence`. It is ordered **after**
  `multi_omics_convergence` and **before** `emerging`. The gate is the
  candidate's own layer is significant AND at least one linked layer is
  significant. High-priority and shared-identity gates are evaluated first, so
  a link never promotes a candidate past a stronger tier and never demotes one.
- **Score:** linked layers do **not** enter `n_supporting_layers` or
  `assay_diversity_score`. The transparent score stays a function of
  shared-identity evidence only. The tier and the card carry the link. (See
  decision B for the alternative.)
- Direction across the link is reported but not required to agree. That is the
  same policy as shared-identity convergence. A product rising while its enzyme
  falls has several valid explanations.

## 6. Evidence card

There is a new subsection under "Multi-omics context":
**"Metabolite–gene links (biochemical, weaker than shared gene identity)"**.
It contains:

- A plain statement: "Linked through KEGG reaction annotation, not measured.
  The enzyme assignment is inferred by KEGG orthology. The link means the gene's
  product is annotated to catalyse a reaction involving this metabolite. It does
  not mean the gene regulates or controls the metabolite's level."
- A table with study, layer, linked feature, compound, reaction, KO, side,
  `link_confidence`, n genes for the compound, direction, adj. p, same
  phenotype, significant, and counts.
- The source line: link table file, KEGG release, access date, sha256.
- Excluded links listed with their reason (`currency_metabolite`,
  `hub_compound`), so the filter is visible rather than silent.
- A stressor caveat when the linked records come from a different stressor
  (for example, a hypoxia metabolite and a heat transcript under
  `growth_under_stress`).

## 7. Decisions (confirmed 2026-10-16)

**A. Source and licensing: builder only, no committed KEGG table.** KEGG
REST is provided for academic use, and KEGG does not clearly grant
redistribution of derived tables. AREE therefore commits the builder
(`aree build-metabolite-links`), the curated compound map, the currency list
and the provenance schema (`schemas/metabolite_gene_links_provenance.schema.json`).
It does not commit the KEGG-derived link table. `data/mappings/metabolite_gene_links_mgigas.*`
is git-ignored. When it has not been built, real-data candidates carry
`link_table_status = not_built` and no linked convergence: the absence is
reported, never silent. CI uses only the simulated demo table. An open-licence
route (Rhea + ChEBI + an eggNOG-mapper run on the RS_2024_06 proteome) remains
future work, because it would allow the table to be committed.

**B. Scoring: tier and card only.** Linked layers never enter
`n_supporting_layers` or `assay_diversity_score`.

**C. A fourth tier**, `metabolite_linked_convergence`, after
`multi_omics_convergence` and before `emerging`.

**D. Demo:** two fabricated demo genes were added to the demo crosswalk and to
GIGAS_HEAT01 `chronic_heat_vs_control`:

- `ldh` (LOC105342001): significant, linked to lactate through R00703/K00016.
- `sdha` (LOC105342002): not significant, linked to succinate through R02164/K00234.

Both are simulated. The KEGG reaction and KO IDs in the demo table are real
KEGG identifiers, but the gene IDs they point to are fabricated.

**E. Cap:** `MAX_GENES_PER_COMPOUND = 25` (`src/mappings/metabolite_links.py`).

### Confidence rule as implemented

A metabolite→gene link is never `exact` or `one_to_one_ortholog`. A compound
that reaches more than one crg gene is `one_to_many_ortholog`. A compound that
reaches exactly one gene is `inferred`, because the enzyme assignment is still
KEGG orthology transfer, not oyster biochemistry. A name→compound match by
`curator_judgement` makes the link `inferred`, which is the weaker rank in
`mapping_confidence.yaml`. The same applies when one reported name maps to
more than one compound. For example, `lactate` maps to both L-lactate
(C00186) and D-lactate (C00256), because LC-MS does not resolve
stereoisomers, and molluscs carry both L- and D-lactate dehydrogenases.

KEGG writes every reaction equation with `<=>`. The table therefore records
the compound's `equation_side` (`left`/`right`), not substrate or product: the
equation does not encode the physiological direction.

Pumps and transporters have no KEGG REACTION, so they get no link. For
example, SERCA (K05853) has none. A significant ATP change and a significant
`atp2a` change are therefore not linked by this mechanism, and ATP is a
currency metabolite in any case.

### First local build (2026-10-16, not committed)

`aree build-metabolite-links` against KEGG REST produced 191 links for 6 of the
7 queried compounds:

- ATP was skipped as a currency compound.
- Pyruvate (47 genes) exceeds the cap and is a hub, so it is not credited.
- 10,875 of 10,879 KEGG `crg` genes with a KO resolved to current RS_2024_06
  GeneIDs. The 4 unresolved IDs are rRNA/tRNA pseudo-IDs (`s-rRNA-1`,
  `tRNA-Arg`, ...), and they are listed in the provenance.
- One biologically meaningful result: **no `crg` gene carries the
  L-lactate dehydrogenase KO (K00016)**, so (S)-lactate has no link, and
  `lactate` links only through D-lactate (C00256, 2 genes). This is consistent
  with bivalve anaerobic metabolism relying on D-LDH and opine dehydrogenases.
  It is also why reported `lactate` is mapped to both stereoisomers rather
  than assumed to be L. The simulated demo, by contrast, links lactate through
  K00016 to a fabricated gene. That is a demo device, not oyster biology.

## 8. Assumptions

- KEGG's crg KO assignments are adequate as *inferred* enzyme annotation for
  *C. gigas*. They come from orthology, not oyster biochemistry.
- Reaction participation is the only link type in the MVP. Pathway
  co-membership is deliberately **not** a link, because it would connect
  almost everything.
- The real table is built for RS_2024_06 only. Other assemblies resolve
  through NCBI GeneIDs, as identifiers already do.
- The only real metabolomics study (DUPERRET2025_POMS_STARVE) failed QC
  (sr320/AREE#70). Real-data convergence is therefore not expected yet, and
  tests use the demo.
