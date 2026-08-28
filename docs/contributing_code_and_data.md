# Contributing code and curated datasets

AREE accepts two distinct kinds of contribution: code changes, and curated
public-dataset additions. They go through the same pull-request mechanism
but have different review concerns.

## Contributing code

For coding conventions, test expectations, linting, and the PR process for
changes to `src/`, `schemas/`, workflow scaffolds, or the app, see
[CONTRIBUTING.md](../CONTRIBUTING.md) at the repository root. This page does
not duplicate that guide — it exists specifically to cover the
curated-dataset contribution path, which is distinct enough from ordinary
code contribution to warrant its own explanation.

## Contributing a curated dataset

Adding a real public study to the registry is the main way AREE grows beyond
its demo data. A dataset contribution PR typically includes:

1. **A study registration YAML** at `registry/studies/STUDY_ID.yaml`, filled
   in per [adding_a_study.md](adding_a_study.md). Do **not** set
   `simulated: true` for a real dataset, and do fill in real `accessions`
   (DOI/BioProject/GEO/SRA/ENA/ProteomeXchange) — see
   `schemas/study.schema.json` for the accepted fields.
2. **A results file** for each comparison declared in the study — a
   processed results table (DE table, DMR table, protein/metabolite
   abundance table) placed under `data/` in a location referenced by
   `comparisons[].results_file`, or, for `raw_reanalysis` studies, evidence
   that the appropriate workflow scaffold in `workflows/` was run and its
   standardized output staged.
3. **Validation output** — run and pass:
   ```bash
   aree validate-study registry/studies/STUDY_ID.yaml
   ```
   before opening the PR, and include the harmonization step if the results
   file is ready:
   ```bash
   aree harmonize --study STUDY_ID
   ```
4. **Honest metadata about data quality and completeness.** If raw QC could
   not be independently verified (processed-results-only mode), if
   replication is low, or if the phenotype definition is ambiguous, set the
   corresponding entries in `quality_flags`
   (`registry/controlled_vocabularies/quality_flags.yaml`) and describe the
   issue in `limitations` rather than omitting it. See
   [governance_and_provenance.md](governance_and_provenance.md) — AREE treats
   missingness and quality caveats as required, reportable information, not
   something to leave implicit.
5. **A PR description** stating: the source citation/accession, which
   `analysis_mode` was used and why, and any identifier-mapping caveats you
   are aware of (e.g. the study only reports legacy gene symbols).

## What a reviewer will check

- Schema validity (`aree validate-study` passes)
- No duplicate `study_id`
- `resilience_classification` on each comparison is defensible given the
  phenotype (see [resilience_vs_exposure.md](resilience_vs_exposure.md)) —
  this is a scientific judgment call, not just a schema check
- `stressor_original` preserves the source study's actual wording rather than
  being pre-collapsed into the standardized term
- `genome_assembly` resolves to an entry in
  `data/reference/genome_assemblies.yaml`; if the assembly is new, add it there
  with an accession you verified yourself against NCBI, not one copied from
  another entry (see [handling_genome_versions.md](handling_genome_versions.md)).
  `aree validate-study` fails until you do
- `species` resolves to a term in
  `registry/controlled_vocabularies/species.yaml`. If validation warns that you
  used an accepted synonym, that is fine — the reported name is preserved — but
  check the term lists the synonym you used
- Provenance fields (`provenance.registered_by`, `provenance.date_registered`)
  are filled in

## Proposing vocabulary or schema changes

If a dataset does not fit the existing phenotype/stressor/tissue/life-stage
vocabularies, propose a new controlled-vocabulary term as part of the same PR
rather than working around the mismatch — see
[defining_a_phenotype.md](defining_a_phenotype.md) for the process. Schema
changes (`schemas/*.json`) are a higher bar and should be raised as an issue
first, since they affect every existing registered study.

## Related documentation

- [adding_a_study.md](adding_a_study.md)
- [resilience_vs_exposure.md](resilience_vs_exposure.md)
- [governance_and_provenance.md](governance_and_provenance.md)
- [CONTRIBUTING.md](../CONTRIBUTING.md) — code contribution guide
