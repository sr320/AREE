# Handling genome versions

Genome assembly and annotation drift is one of the main reasons public
omics results become hard to compare over time. AREE tracks both explicitly
rather than assuming "the oyster genome" is a single stable thing.

## The fields

Every study record (`schemas/study.schema.json`) carries:

- `genome_assembly` (required) — the reference genome assembly identifier
  used for alignment/annotation in that study, e.g. `"cgigas_uk_roslin_v1"`.
- `annotation_version` (optional but strongly recommended) — the specific
  gene annotation build, since the same assembly can have multiple
  annotation releases over time, e.g. `"cgigas_annotation_v1_demo"` or, for a
  deliberately mismatched demo case,
  `"legacy_oyster_v9_symbols_only_demo"` (`registry/studies/GIGAS_LARV05.yaml`
  — a study annotated with an older, symbol-only gene set).

Every evidence record (`schemas/evidence.schema.json`) independently carries
its own `genome_assembly` and `annotation_version`, copied from the study/
comparison that produced it, so a downstream consumer of the evidence table
never has to join back to the study registry just to know which reference
was used for a given effect estimate.

## Why this is not collapsed into one global "current genome"

Because studies span years, a resilience evidence engine that only recorded
"the gene" without recording which assembly/annotation defined that gene
would silently misattribute results whenever locus definitions shift between
annotation releases. Keeping `genome_assembly` + `annotation_version` on
every record means:

- identifier harmonization (`src/harmonize/identifiers.py`) can be extended
  to select the correct crosswalk for the assembly/annotation actually used,
  rather than assuming one crosswalk fits everything (see
  [adding_a_species.md](adding_a_species.md) for the current single-crosswalk
  limitation), and
- a reviewer can immediately see if two "matching" identifiers across studies
  were actually called against different annotation versions, which is a
  legitimate reason to treat a match with more caution even at
  `mapping_confidence: exact`.

## `data/reference/genome_assemblies.yaml`

This file documents assembly metadata referenced by studies:

```yaml
assemblies:
  - assembly_id: "cgigas_uk_roslin_v1"
    species: "Crassostrea gigas"
    ncbi_assembly_accession: "PLACEHOLDER_CONFIRM_BEFORE_REAL_USE"
    notes: "Chromosome-level assembly label used as the harmonization target for all demo studies. Demo DMR coordinates (data/demo/methylation/) are illustrative and not derived from real alignment or a verified real assembly accession."
```

**This is explicitly a placeholder, not a verified real accession.** The
file's own header comment states: "`assembly_id` and `ncbi_assembly_accession`
below are PLACEHOLDERS for demo purposes only. Do not treat
`ncbi_assembly_accession` as a verified real accession — when curating a real
study, look up and confirm the actual current assembly accession from
NCBI/Ensembl yourself rather than reusing this value." Do not cite this value
in any real study registration or publication; it exists only so the demo
has a self-consistent label to point at.

## Adding a real assembly entry

When curating a real (non-demo) study:

1. Look up the actual assembly accession from NCBI or Ensembl yourself (do
   not copy the demo placeholder).
2. Add a new entry to `data/reference/genome_assemblies.yaml` with a
   descriptive `assembly_id`, the confirmed `ncbi_assembly_accession`, and
   `notes` on annotation provenance.
3. Use that `assembly_id` value consistently in `genome_assembly` across the
   study record and every comparison's harmonized evidence records.
4. If the source study used an older annotation than your reference
   crosswalk, record that in `annotation_version` and expect a lower
   `mapping_confidence` — see [identifier_mapping.md](identifier_mapping.md)
   for how `GIGAS_LARV05`'s legacy symbol-only annotation is handled.

## Handling a genome version change for an already-registered study

If a study needs to be re-harmonized against a newer assembly/annotation
(e.g. the reference was updated after the study was first registered):

1. Update `genome_assembly`/`annotation_version` on the affected study
   record, and log the change in `provenance.curation_notes` (do not edit
   history silently — see [governance_and_provenance.md](governance_and_provenance.md)).
2. Re-run `aree harmonize --study STUDY_ID` to regenerate evidence records
   against the current identifier crosswalk.
3. Because `workflow_version` and `date_generated` are recorded on every
   evidence record, both the old and new harmonization runs remain
   traceable — nothing is destructively overwritten in a way that hides
   which version produced which row.

## Related documentation

- [adding_a_species.md](adding_a_species.md)
- [identifier_mapping.md](identifier_mapping.md)
- [governance_and_provenance.md](governance_and_provenance.md)
