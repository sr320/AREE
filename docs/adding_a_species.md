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
   is `is_ncbi_reference: true`; a test enforces it. The exception is a species
   with no annotated assembly yet (the basket cockle): its assembly is
   registered with `is_ncbi_reference: false` and a null crosswalk, even when
   NCBI calls it the reference genome. Validation cross-checks the
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

## Species identity and reads mapped to another genome

Two rules protect the species label itself:
- **Hybridizing complex.** A species whose term has a `complex` (the *Mytilus edulis* complex) needs
  `species_identity` on every study, and any basis short of `genotyped` carries
  `species_assignment_unverified`.
- **Reference-species mismatch.** A study whose reads were mapped to another species' genome must
  declare `reference_species_mismatch`. Its identifiers follow the assembly, but its evidence is
  partitioned under the study animal's taxid.

Both are explained, with the validation rules and how the dashboard shows them, in
[species_identity.md](species_identity.md).

## Worked example: the Manila clam

The integration plan proposed the first mussel as the worked example. No real
mussel study has been harmonized yet (the only mussel study is the simulated
`MYTGAL_HEAT01`), so this follows the Manila clam, *Ruditapes philippinarum*,
the first species other than the oyster to go through every step on real
data. Dates are when each step happened.

1. **Vocabulary term (2026-10-11).** `ruditapes_philippinarum`, taxid 129788,
   with seven accepted synonyms from NCBI Taxonomy, among them *Venerupis
   philippinarum* and *Tapes philippinarum*. Studies use all of them, and a
   missing synonym would have split the species silently.
2. **Assembly (2026-10-11).** `ASM2657151v2` (GCF_026571515.1, annotation
   GCF_026571515.1-RS_2023_10, University of Bologna). It is the only
   RefSeq-annotated Manila clam assembly, and it is contig-level (N50 about
   183 kb). Its `notes` say so, because gene-level evidence is unaffected but
   region-level evidence would be.
3. **Crosswalk (2026-10-11).**

   ```bash
   aree build-crosswalk --taxid 129788 --slug rphil \
     --organism "Ruditapes philippinarum, Manila clam"
   ```

   This gave 48,038 genes, 57.1% with an Ensembl xref. Only 18 genes link to UniProt, because
   UniProt's clam entries come from GenBank submissions, not the RefSeq annotation.
4. **Reference bundle and proteome (2026-10-12).**

   ```bash
   aree fetch-reference --assembly ASM2657151v2
   aree fetch-reference --assembly ASM2657151v2 --proteins
   ```

   The bundle held 65,336 transcripts and 47,388 genes, every one of them in the crosswalk. The
   proteome went into the first orthogroup run (`orthofinder3.1.5_2026-10-13`), which shares
   85-88% of each clam's genes with at least one other species.
5. **First study (2026-10-12 to 13).** `aree start-study --accession
   PRJNA1132140 --study-id TUCCI2025_HEAT_PRIMING`, then `aree
   fetch-samplesheet ... --library-strategy RNA-Seq`. Three things the
   deposit did not show came from the paper:
   - **Double-counted groups.** The BioProject also held 16S amplicon runs from the same animals,
     which doubled every group until `--library-strategy` existed.
   - **3′ Tag-seq libraries.** The libraries were QuantSeq 3′ Tag-seq although ENA said PolyA. That
     needed the workflow's `--tag_seq` mode, with no length correction anywhere.
   - **Tanks, not clams, as the unit.** Tank was the unit of treatment, with tank IDs not deposited.
     The study was fitted on clams with `tank_level_replication`, an explicit decision recorded in
     the YAML.
6. **Second study (2026-10-15 to 17).** `PERUZZA2023_HEATWAVE` (PRJNA906528):
   - **Library type.** The same library type, confirmed from the paper and from a one-sided poly-A
     signal in the reads.
   - **Independence.** Confirmed from purchase dates: clams bought a year apart from the same
     hatchery.
   - **Sex.** Fitted as a covariate because the arms were unbalanced.
7. **First pool.** Both studies' heatwave-versus-control comparisons are
   `meta_analysis_primary` on `thermal_tolerance`, classified
   `exposure_only`. That made the first within-species clam pool, k = 2.
8. **Across species.** With two oyster thermal-tolerance studies already on
   `main`, oyster and clam then pooled through 5,975 single-copy
   orthogroups. The dashboard gate opened for thermal tolerance. The
   replication gate was met in 5,633 orthogroups, but no row converges:
   both clam comparisons are exposure-only.

What the example teaches:

- **Read the methods before the metadata.** Library type, replication unit
  and survivor sampling were all invisible in the deposit.
- **Prespecify the pooled comparison before running**, and record it in
  the YAML with the date. Both clam studies did.
- **A species' first pool can be exposure-only.** The pipeline will say so
  on every card, and cross-species convergence on a tolerance phenotype
  then needs a study that measured tolerance.

### A species NCBI has not annotated: *Mytilus coruscus*

*M. coruscus* (taxid 42192) has no NCBI Gene records. Its only annotated
assembly, MCOR1.1 (GCA_011752425.2), carries a GenBank annotation from its
authors. It was registered with `is_ncbi_reference: false` and its crosswalk
built by the GFF route (step 3 above): 58,540 genes anchored on the authors'
locus tags, and 9,877 with a Swiss-Prot hit from a recorded DIAMOND search. It
is not in an orthogroup run yet, and no study is registered.

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
