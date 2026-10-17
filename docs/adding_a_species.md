# Adding a species

AREE began as a single-species (*Magallana gigas*) resource, but `species` is a
plain data field on every study and evidence record, not an assumption baked
into the schemas or code. Onboarding another species is a data-addition
exercise. The mussel and clam work in
[mussel_clam_integration_plan.md](mussel_clam_integration_plan.md) is the
running example: six mussel and clam species, seven of their assemblies, and a
crosswalk for each species are registered, so their studies can be harmonized.

## What is already species-agnostic

- `schemas/study.schema.json` and `schemas/evidence.schema.json` both carry
  `species`, validated against `registry/controlled_vocabularies/species.yaml`.
- The phenotype and stressor ontologies are defined in organism-agnostic terms
  (survival, thermal tolerance, temperature, salinity, hypoxia, and so on).
- Meta-analysis groups by `feature_id_standardized` + `phenotype` +
  `feature_type` + `simulated` + **`species_taxid`**. The species key keeps two
  species from pooling into one estimate, so same-feature evidence from two
  species produces two pooled rows, by design. Combining across species is only
  meant to happen through `orthogroup_id`, which harmonization now attaches from
  the active orthogroup run (see
  [identifier_mapping.md](identifier_mapping.md#orthogroups)).
- **Identifier resolution follows each study's assembly.** Harmonizing a real
  study reads the `identifier_crosswalk` registered for its `genome_assembly`
  in `data/reference/genome_assemblies.yaml`. Studies of different species in
  one run (`make real-pool`, the Pages build) each resolve against their own
  crosswalk, with no environment variable or CI edit. See
  [identifier_mapping.md](identifier_mapping.md#which-crosswalk-is-in-force).

## What you actually need to add

1. **A species vocabulary term.** Add a term to
   `registry/controlled_vocabularies/species.yaml` with the accepted
   `scientific_name`, its `ncbi_taxid`, and any `accepted_synonyms`.

   Take the synonyms seriously. The Pacific oyster is the cautionary example:
   it was moved from *Crassostrea* to *Magallana*, both names are in current
   use, and both resolve to taxid 29159 so they group as one animal. Omit a
   synonym and you get a species silently split in two, with no error, just
   two half-powered meta-analyses. The Manila clam has the same problem with
   *Ruditapes*, *Venerupis* and *Tapes*.

   Add only synonyms that NCBI Taxonomy lists for the taxid. E-utilities
   returns them under `OtherNames`:

   ```bash
   curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=129788"
   ```

   If the species belongs to a group that hybridizes and is often mislabelled
   in public deposits, set `complex` on every member (the three *Mytilus
   edulis*-complex mussels carry `complex: mytilus_edulis_complex`). See
   *Species in a hybridizing complex* below for what that requires of studies.

   `aree validate-study` rejects a species that is not in this file, and warns
   (without failing) when a study uses a synonym rather than the accepted name.

2. **A `genome_assembly` reference.** Add an entry to
   `data/reference/genome_assemblies.yaml` for each assembly studies of the
   species use, with `assembly_id`, verified RefSeq and GenBank accessions,
   `ncbi_taxid`, annotation release, the date you `verified` it, and
   `identifier_crosswalk` (null until step 3). Exactly one assembly per species
   is `is_ncbi_reference: true`; a test enforces it. Validation cross-checks the
   assembly taxid against the species. See
   [handling_genome_versions.md](handling_genome_versions.md); look accessions
   up yourself against the NCBI Datasets API.

3. **An identifier crosswalk.** For a species with an NCBI RefSeq annotation
   (a `GCF_` accession with an annotation release), build one from NCBI Gene and
   UniProtKB:

   ```bash
   aree build-crosswalk --taxid 129788 --slug rphil \
     --organism "Ruditapes philippinarum, Manila clam"
   ```

   This writes `data/reference/crosswalk/<slug>_gene_id_crosswalk.tsv`, a
   retired-GeneID sidecar, and a provenance JSON recording the taxid. Commit
   all three, then set `identifier_crosswalk` to the TSV path on every assembly
   of that species in `genome_assemblies.yaml`. A test checks that each
   registered crosswalk's provenance taxid matches its assembly's.

   Until this step is done, real studies of the species can be **registered**
   (validation warns) but not **harmonized** (harmonization stops with the
   command above in its error message).

   A species with only a GenBank assembly and an author-supplied annotation has
   no NCBI Gene records. Build its crosswalk from the authors' GFF3 and
   proteins with a recorded DIAMOND search against one UniProtKB release:

   ```bash
   aree build-crosswalk --taxid 42192 --slug mcoru \
     --from-gff genomic.gff --proteins protein.faa \
     --uniprot-fasta uniprot_sprot.fasta.gz \
     --uniprot-release "UniProtKB/Swiss-Prot Release 2026_03 of 02-Sep-2026" \
     --uniprot-url https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz \
     --diamond /path/to/diamond --threads 6
   ```

   Register the assembly first with `is_ncbi_reference: false` and its
   `identifier_crosswalk` set; the builder reads the annotation description
   from it. Such a species has no NCBI reference assembly, and tests allow
   that only for a GFF-route crosswalk. Commit the crosswalk and its
   provenance; there is no retired-GeneID sidecar. See
   [identifier_mapping.md](identifier_mapping.md#genomes-without-ncbi-gene-records-gff-route-wp7).

4. **Species-appropriate tissue, life-stage, phenotype and stressor terms**,
   only if the existing vocabulary does not fit. Add new terms rather than
   repurposing existing ones. The mussel and clam additions are an example:
   `hemocytes`, `foot`, `byssal_gland` and `siphon` tissues; the
   `emersion_desiccation` stressor; and the `byssal_attachment` and
   `burrowing_performance` phenotypes, each with a `resilience_relevance`.

5. **Orthology context across species.** Fetch the species' proteins
   (`aree fetch-reference --assembly <acc> --proteins`), then rebuild the
   orthogroup run so it includes the new species (`aree build-orthogroups`)
   and make it active. Orthogroups are inferred jointly, so adding a species
   means a new run, not an addition to the old one. See
   [identifier_mapping.md](identifier_mapping.md#orthogroups).

## Species in a hybridizing complex

Mussels of the *Mytilus edulis* complex (*M. edulis*, *M. galloprovincialis*,
*M. trossulus*) hybridize where their ranges meet, and public deposits are
labelled by the depositor's assumption rather than by genotype. For any species
whose term has a `complex`, `aree validate-study` requires the study to say how
its animals were assigned:

```yaml
species: "Mytilus galloprovincialis"
species_identity:
  basis: collection_site   # genotyped | morphology | collection_site | authors_statement | unstated
  notes: "Adriatic longline farm stock; no species-diagnostic markers reported."
quality_flags:
  - species_assignment_unverified
```

Any basis other than `genotyped` must carry the quality flag
`species_assignment_unverified`, which then travels with every evidence record.
Do not drop the flag to make a pool larger. The complex is a display and
curation grouping only: evidence is still partitioned by `species_taxid`.

## Reads mapped to another species' genome

Some species have no assembly of their own (*Mytilus chilensis*, the greenshell
mussel), and studies map their reads to a relative's genome. Validation fails
when a study's species and its assembly are different taxa, unless the study
declares the quality flag `reference_species_mismatch` and explains it in
`limitations`. The evidence is then partitioned under the study animal's
`species_taxid`, while identifiers resolve against the assembly's crosswalk,
because the gene IDs really are the reference species'.

## What you do not need to change

- No schema changes: `species`, `species_identity` and `species_taxid` are
  already in place.
- No changes to meta-analysis, scoring, harmonization code, the Makefile, or CI
  workflows. Crosswalk selection is driven by `genome_assemblies.yaml`.

## Registering the first study for a new species

Follow [adding_a_study.md](adding_a_study.md) as normal, setting `species` and
`genome_assembly` to the new values. `aree validate-study` checks the species
and assembly exist and agree, enforces the species-complex and
reference-mismatch rules above, and warns if the assembly has no crosswalk yet.

## Related documentation

- [mussel_clam_integration_plan.md](mussel_clam_integration_plan.md)
- [handling_genome_versions.md](handling_genome_versions.md)
- [identifier_mapping.md](identifier_mapping.md)
- [roadmap.md](roadmap.md)
- [design.md](design.md#9-explicit-assumptions-mvp-scope)
