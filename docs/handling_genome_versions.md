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

This file is the registry of assemblies AREE knows about, and it is
**load-bearing**: `aree validate-study` rejects a `genome_assembly` value that
does not resolve to an entry here, so a typo or an undocumented assembly is
caught at registration rather than surfacing later as evidence that cannot be
compared. It also cross-checks the assembly's taxid against the study's species.

Every accession in the file was verified against the NCBI Datasets API on the
date recorded in its `verified_against` header. Two assemblies are registered
for the Pacific oyster:

| assembly_id | accession | role |
|---|---|---|
| `xbMagGiga1.1` | `GCF_963853765.1` | current NCBI reference (Sanger); the annotation AREE's crosswalk is built from |
| `cgigas_uk_roslin_v1` | `GCF_902806645.1` / `GCA_902806645.1` | Roslin assembly; what most pre-2024 studies, including `HESSER2024_VCOR`, were mapped against |

A study may write the `assembly_id`, either accession, or the two together as
publications often do (`"GCA_902806645.1 (cgigas_uk_roslin_v1)"`) — all resolve
to the same record.

## The crossing between a study's assembly and AREE's annotation

This is the case the file exists to make visible, and it is the normal case,
not an edge case.

`HESSER2024_VCOR` was mapped by its authors against the Roslin assembly. AREE
standardizes its identifiers against the **current NCBI annotation**
(`xbMagGiga1.1`, release `GCF_963853765.1-RS_2024_06`), because that is what
the crosswalk is built from. Those are different assemblies from different
sequencing centres.

That is not an error. NCBI Gene IDs are stable identities that persist across
assemblies and re-annotations, which is exactly why the identifier hierarchy
prefers them — see [identifier_mapping.md](identifier_mapping.md). But a reader
must be able to see that the crossing happened, so every evidence record now
carries both:

- `genome_assembly` — the assembly the **source study** used;
- `identifier_annotation_release` — the annotation `feature_id_standardized`
  actually refers to.

For the real study those read `GCA_902806645.1 (cgigas_uk_roslin_v1)` and
`GCF_963853765.1 (xbMagGiga1.1), annotation release GCF_963853765.1-RS_2024_06`
respectively. For demo studies `identifier_annotation_release` is **empty**: the
demo crosswalk is synthetic and has no annotation provenance to report, and
inventing one would be worse than leaving it blank.

When the two disagree, treat coordinate-based evidence (methylation regions,
QTL intervals) with far more caution than gene-level evidence: coordinates do
not survive an assembly change, whereas Gene IDs largely do. AREE does not yet
perform liftover — see [roadmap.md](roadmap.md).

## Adding a real assembly entry

When curating a real (non-demo) study:

1. Look up the actual assembly accession from NCBI or Ensembl yourself. The
   NCBI Datasets API answers this directly:

   ```bash
   curl -s "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_963853765.1/dataset_report"
   ```

2. Add an entry to `data/reference/genome_assemblies.yaml` with a descriptive
   `assembly_id`, the confirmed accessions, `ncbi_taxid`, and `notes` on
   annotation provenance. Registration will fail until you do — that is
   deliberate.
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
