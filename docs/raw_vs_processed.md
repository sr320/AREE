# Raw reanalysis vs. processed-results harmonization

Not every public dataset comes with usable raw data. AREE supports two modes
so that a study without deposited FASTQ/BAM/spectra can still contribute
evidence, while being explicit about the difference in what could be
independently verified.

## The `analysis_mode` field

Every study declares exactly one value in `analysis_mode`
(`schemas/study.schema.json`):

- **`raw_reanalysis`** — raw data (FASTQ, BAM, WGBS reads, raw spectra, etc.)
  is available and is (or will be) run through the standardized Nextflow
  workflow scaffolds in `workflows/` (`rnaseq/`, `methylation/`,
  `proteomics/`, `metabolomics/`), producing the standardized result tables
  that `aree harmonize` then converts into evidence records.
- **`processed_results_harmonization`** — only a processed/summary results
  table is available (a DE table, a DMR table, a protein abundance table, a
  metabolite feature table). `src/harmonize` maps that table directly into
  the same evidence schema without re-running upstream QC/alignment/calling,
  because there is nothing upstream to re-run.

## How to pick

Use `raw_reanalysis` only when raw data is actually accessible and you intend
to (or already have) run it through a workflow scaffold end to end. Use
`processed_results_harmonization` whenever:

- the source deposits only a supplementary results table (common for older
  studies or proteomics/metabolomics submissions),
- raw data exists but is restricted/embargoed (`data_availability.status`),
- or you have not yet run the raw-data pipeline and want the study registered
  and its existing published results harmonized in the meantime.

In the demo set, `GIGAS_HEAT01`, `GIGAS_OA02`, and `GIGAS_LARV05` (all
RNA-seq) and `GIGAS_PATH03` (methylation) use `raw_reanalysis`, reflecting
that their standardized result tables were produced end-to-end for the demo.
`GIGAS_SAL04` (proteomics) and `GIGAS_GROW06` (metabolomics) use
`processed_results_harmonization` — both note in `data_availability.access_notes`
that only a processed table was ever available. `GIGAS_SAL04`'s
`quality_flags` explicitly include `processed_only` for this reason.

## What changes, and what does not

What differs:

- Raw mode is expected to produce (or has produced) the workflow-level
  outputs described in [architecture.md](architecture.md#layer-2) — QC
  metrics, alignment/calling logs, a machine-readable manifest with tool
  versions and checksums.
- Processed mode cannot independently verify upstream QC. `aree harmonize`
  does not fabricate QC metrics it cannot compute; instead, studies
  registered this way should carry the `processed_only` quality flag
  (`registry/controlled_vocabularies/quality_flags.yaml`) so downstream
  consumers of the evidence table know raw QC was unverifiable.

What does **not** differ: both modes converge on the exact same evidence
schema (`schemas/evidence.schema.json`). Meta-analysis and prioritization
code (`src/meta_analysis`, `src/prioritize`) read evidence records without
knowing or caring which mode produced them — a `raw_reanalysis` row and a
`processed_results_harmonization` row are pooled together as long as they
share `feature_id_standardized` + `phenotype` + `feature_type`. The only
place the distinction shows up downstream is via `quality_flags`, which do
feed into `quality_score` in candidate scoring (see
[interpreting_candidate_scores.md](interpreting_candidate_scores.md)).

## Running harmonization

Regardless of mode, harmonization is invoked the same way:

```bash
aree harmonize --study STUDY_ID --input path/to/results.tsv
```

or, to harmonize every comparison already associated with a study's
`results_file` entries:

```bash
aree harmonize --study STUDY_ID
```

Output is appended to the shared evidence table at
`reports/evidence/evidence_table.tsv`.

## Related documentation

- [architecture.md](architecture.md)
- [adding_a_study.md](adding_a_study.md)
- [design.md](design.md#7-raw-vs-processed-results-modes)
- [roadmap.md](roadmap.md) for the current status of raw-mode Nextflow execution
