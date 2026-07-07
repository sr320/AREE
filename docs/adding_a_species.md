# Adding a species

AREE's MVP demo is single-species (*Crassostrea gigas*), but `species` is a
plain data field on every study and evidence record — not a hardcoded
assumption baked into the schemas or code. Onboarding a second species is a
data-addition exercise, not a schema change.

## What is already species-agnostic

- `schemas/study.schema.json` and `schemas/evidence.schema.json` both carry
  `species` as a free-text field (e.g. `"Crassostrea gigas"`); nothing in
  the schema restricts it to oysters.
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
- Meta-analysis and prioritization code (`src/meta_analysis`, `src/prioritize`)
  group and score by `feature_id_standardized` + `phenotype` + `feature_type`
  and do not reference species directly, so nothing there needs to change.

## What you actually need to add

1. **A `genome_assembly` reference.** Add an entry to
   `data/reference/genome_assemblies.yaml` for the new species' assembly,
   following the existing structure (`assembly_id`, `species`,
   `ncbi_assembly_accession`, `notes`). See
   [handling_genome_versions.md](handling_genome_versions.md) — verify the
   real accession yourself; do not copy the current placeholder value.
2. **A new identifier crosswalk**, most likely
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
3. **Species-appropriate tissue/life-stage terms**, if the existing bivalve
   vocabulary doesn't fit (e.g. a finfish species will need different tissue
   terms than `gill`/`mantle`/`digestive_gland`). Add new terms to
   `tissue_types.yaml` / `life_stages.yaml` rather than repurposing existing
   ones.
4. **Orthology context across species**, if you want cross-species candidate
   comparison (e.g. linking a *C. gigas* candidate to an ortholog in a second
   shellfish species). The `orthogroup_id` field exists in the crosswalk and
   evidence schema for exactly this purpose, but the demo crosswalk's
   orthogroup IDs are illustrative single-species placeholders, not the
   output of a real cross-species orthology caller (see
   [roadmap.md](roadmap.md) — OrthoFinder-based mapping is future work).

## What you do not need to change

- No schema changes: `schemas/study.schema.json` and
  `schemas/evidence.schema.json` already carry `species` as a required
  field.
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
