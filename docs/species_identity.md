# Species identity and reference-species mismatch

AREE keeps every species' evidence apart: pooling, ranking and the
false-discovery families are all partitioned by `species_taxid`. Species are
compared only through orthogroups, and only behind the dashboard gate
([interpreting_meta_analysis.md](interpreting_meta_analysis.md#cross-species-results---level-orthogroup)).
That separation is only as good as the species label on each study. Two
situations weaken it, and AREE makes both explicit instead of guessing.

## When the label itself is uncertain: hybridizing complexes

The blue mussels *Mytilus edulis*, *M. galloprovincialis* and *M. trossulus*
hybridize wherever their ranges meet. In much of the North Atlantic, the
Baltic, the Pacific Northwest and the Southern Hemisphere, a mussel's species
cannot be told from its shell. Public deposits are labelled by the depositor's
assumption: often the species usually farmed in that region, rarely a
genotype. A study labelled *M. edulis* may be partly or wholly *M. trossulus*,
or hybrids.

The vocabulary marks these species with
`complex: mytilus_edulis_complex` (`registry/controlled_vocabularies/species.yaml`).
For any species with a `complex`, `aree validate-study` requires the study to
say how its animals were assigned:

```yaml
species: "Mytilus edulis"
species_identity:
  basis: collection_site   # genotyped | morphology | collection_site | authors_statement | unstated
  notes: "Outer Oslofjord; no species-diagnostic markers reported."
quality_flags:
  - species_assignment_unverified
```

| `basis` | Meaning | Needs `species_assignment_unverified` |
|---|---|---|
| `genotyped` | Species-diagnostic markers or sequence data on the study animals themselves | no |
| `morphology` | Shell or anatomical characters | yes |
| `collection_site` | Inferred from where the animals were collected | yes |
| `authors_statement` | The paper or deposit names the species with no stated basis | yes |
| `unstated` | Nothing recorded | yes |

Validation fails if a complex species has no `species_identity`, or if a basis
short of `genotyped` lacks the flag. The flag then travels with every evidence
record the study produces. It counts against the candidate's quality score and
is listed under Limitations on its evidence cards.

What this does and does not do:

- **It does not relabel anything.** Evidence stays under the taxid the study
  declared. AREE has no way to re-assign species after the fact, and guessing
  would be worse than a flagged label.
- **It does not merge the complex.** `complex` is a display and curation
  grouping. *M. edulis* and *M. galloprovincialis* evidence never pool at gene
  level, exactly as two unrelated species would not. Cross-species pooling
  through orthogroups treats them as different species too.
- **It is visible.** The dashboard marks complex species with † in the species
  table and adds a *Species identity* column to the studies table whenever a
  listed study belongs to a complex. The Streamlit app shows
  `species_identity_basis` in the study registry and the complex in its
  species table.

When curating: do not drop the flag to enlarge a pool. If the paper genotyped
its animals (mitochondrial or nuclear species-diagnostic markers, or SNP
panels), record `genotyped` and name the markers in `notes`. If it reports a
mixture, the honest registration is the mixture's dominant species with
`basis: genotyped`, the flag, and the proportions in `notes` and
`limitations`.

## When the reads were mapped to another species' genome

Some species have no assembly of their own (*Mytilus chilensis*, for
example), and their studies map reads to a relative's genome. The evidence is then about one animal but
its gene identifiers belong to another species.

AREE handles this explicitly:

- **Validation fails by default** when a study's species and its
  `genome_assembly` are different taxa. The study must declare the quality
  flag `reference_species_mismatch` and explain it in `limitations`. Declaring
  the flag when the taxa match is also an error.
- **Identifiers follow the assembly.** The crosswalk is chosen from the
  assembly (`identifier_crosswalk` in `data/reference/genome_assemblies.yaml`),
  because the gene IDs really are the reference species'. See
  [identifier_mapping.md](identifier_mapping.md#which-crosswalk-is-in-force).
- **Evidence is partitioned under the study animal's taxid**, not the
  reference's, so it never pools with the reference species' own evidence at
  gene level.
- **The flag travels with every record**, so a reader sees that a "gene" in
  this study is a read-mapping proxy: divergent genes map poorly, and
  lineage-specific genes have no counterpart at all.

## How to check a registration

```bash
aree validate-study registry/studies/<STUDY_ID>.yaml
```

The errors name the rule and the fix. On the dashboard, the species table
lists every registered species with its crosswalk route and study counts, and
the studies table filters by species.

## Related documentation

- [adding_a_species.md](adding_a_species.md)
- [identifier_mapping.md](identifier_mapping.md)
- [handling_genome_versions.md](handling_genome_versions.md)
- [mussel_clam_integration_plan.md](mussel_clam_integration_plan.md) (WP1, WP8)
