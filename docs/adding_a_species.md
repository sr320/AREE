# Adding a species

AREE's MVP demo is single-species (*Crassostrea gigas*), but `species` is a
plain data field on every study and evidence record — not a hardcoded
assumption baked into the schemas or code. Onboarding a second species is a
data-addition exercise, not a schema change.

## What is already species-agnostic

- `schemas/study.schema.json` and `schemas/evidence.schema.json` both carry
  `species`; nothing in the schema restricts it to oysters. The value is
  validated against `registry/controlled_vocabularies/species.yaml`, so a new
  species needs a term there — see step 1 below.
- The phenotype ontology (`registry/controlled_vocabularies/phenotype_ontology.yaml`)
  and stressor ontology (`registry/controlled_vocabularies/stressor_ontology.yaml`)
  are defined in organism-agnostic terms (survival, thermal tolerance,
  disease resistance, temperature, salinity, hypoxia, etc.) and apply to any
  aquaculture species without modification.
- `tissue_types.yaml` and `life_stages.yaml` are largely bivalve-appropriate
  today (gill, mantle, digestive gland; gamete, embryo, larva, spat) but are
  just YAML lists — add species-appropriate terms as needed rather than
  overloading an existing term to mean something different for a new
  species.
- Meta-analysis groups by `feature_id_standardized` + `phenotype` +
  `feature_type` + `simulated` + **`species_taxid`**. The species key is what
  keeps two species from pooling into one estimate the moment a second species
  is registered, so nothing in `src/meta_analysis` or `src/prioritize` needs to
  change to onboard one — but understand that same-feature evidence from two
  species will produce two separate pooled rows, by design. Combining across
  species is only ever meant to happen through `orthogroup_id`, which is not
  yet populated (see [roadmap.md](roadmap.md)).

## What you actually need to add

1. **A species vocabulary term.** Add a term to
   `registry/controlled_vocabularies/species.yaml` with the accepted
   `scientific_name`, its `ncbi_taxid`, and any `accepted_synonyms`.

   Take the synonyms seriously. The Pacific oyster is the cautionary example
   already in the file: it was moved from *Crassostrea* to *Magallana*, both
   names are in current use, and every study registered in AREE so far says
   *Crassostrea* while the reference crosswalk is built under *Magallana*.
   Because both resolve to taxid 29159, they group as one animal. Omit a
   synonym and you get a species silently split in two — with no error, just
   two half-powered meta-analyses.

   `aree validate-study` rejects a species that is not in this file, and warns
   (without failing) when a study uses a synonym rather than the accepted name.

2. **A `genome_assembly` reference.** Add an entry to
   `data/reference/genome_assemblies.yaml` for the new species' assembly with
   `assembly_id`, verified accessions, and `ncbi_taxid`. Validation cross-checks
   that taxid against the species, so a mismatched pair fails at registration.
   See [handling_genome_versions.md](handling_genome_versions.md) — verify the
   accession yourself against NCBI.
3. **A new identifier crosswalk**, most likely
   `data/mappings/<species>_gene_id_crosswalk.tsv`, following the same column
   structure as `data/mappings/gene_id_crosswalk.tsv`
   (`ncbi_gene_id`, `ensembl_gene_id`, `uniprot_accession`, `locus_id`,
   `gene_symbol`, `orthogroup_id`). `src/harmonize/identifiers.py` currently
   points at a single crosswalk path (`CROSSWALK_PATH`); supporting multiple
   species' crosswalks in the same run requires a small code change to select
   the crosswalk by species (or genome_assembly) rather than a fixed path —
   this is a real, not-yet-implemented extension point, not something you can
   do purely by adding a data file. Track it as implementation work before
   registering a mixed-species batch.
4. **Species-appropriate tissue/life-stage terms**, if the existing bivalve
   vocabulary doesn't fit (e.g. a finfish species will need different tissue
   terms than `gill`/`mantle`/`digestive_gland`). Add new terms to
   `tissue_types.yaml` / `life_stages.yaml` rather than repurposing existing
   ones.
5. **Orthology context across species**, if you want cross-species candidate
   comparison (e.g. linking a *C. gigas* candidate to an ortholog in a second
   shellfish species). The `orthogroup_id` field exists in the crosswalk and
   evidence schema for exactly this purpose, but the demo crosswalk's
   orthogroup IDs are illustrative single-species placeholders, not the
   output of a real cross-species orthology caller (see
   [roadmap.md](roadmap.md) — OrthoFinder-based mapping is future work).

## What you do not need to change

- No schema changes: `schemas/study.schema.json` and
  `schemas/evidence.schema.json` already carry `species` as a required field,
  and `species_taxid` is derived automatically from the vocabulary term.
- No changes to `phenotype_ontology.yaml` or `stressor_ontology.yaml` unless
  the new species genuinely needs a phenotype/stressor concept not already
  covered — these vocabularies were designed to be organism-agnostic (see
  [design.md](design.md#3-controlled-vocabularies)).
- No changes to meta-analysis or scoring code — both operate on standardized
  identifiers and ontology terms, not species names.

## Registering the first study for a new species

Follow [adding_a_study.md](adding_a_study.md) as normal, setting `species`
and `genome_assembly` to the new values. `aree validate-study` will validate
phenotype/stressor/tissue/life-stage terms the same way regardless of
species; it does not currently cross-check that `genome_assembly` matches an
entry in `data/reference/genome_assemblies.yaml`, so keep that file in sync
manually.

## Related documentation

- [handling_genome_versions.md](handling_genome_versions.md)
- [identifier_mapping.md](identifier_mapping.md)
- [roadmap.md](roadmap.md)
- [design.md](design.md#9-explicit-assumptions-mvp-scope)
