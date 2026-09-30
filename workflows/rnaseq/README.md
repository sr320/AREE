# AREE :: RNA-seq workflow

Implements CLAUDE.md Layer 2, section A (RNA-seq): FASTQ QC, adapter
trimming, pseudoalignment/quantification, sample QC, differential
expression, effect-size calculation, and a standardized, enrichment-ready
result table plus a human-readable report.

**Status: executed end to end against real public FASTQ.** On 2026-08-28 this
workflow ran the full `raw_reanalysis` path — FASTQC → fastp → Salmon index →
Salmon quant → MultiQC → DESeq2 → standardize → manifest → report, 38 processes
green — on 11 libraries from BioProject PRJNA1329250 (`CALLA2026_OSHV`), and
produced 22,301 genes with real `lfcSE`. `processed_results_harmonization` mode
also runs end to end.

Two important qualifications:

* the run used **subsampled reads** (first 1M pairs per library) and covered one
  of six comparisons, so it validates the pipeline, not the biology, and its
  output has not been harmonized into the evidence table;
* it ran **natively** via `-profile local` with Homebrew-installed tools. **No
  container has ever been pulled**, so every image tag below remains unverified.

Getting there required fixing sixteen defects, several of which prevented this
workflow from compiling at all on Nextflow 26.04. They are catalogued in
`../../docs/first_raw_reanalysis.md`.

## Two modes

Selected via `params.mode` (default: `processed_results_harmonization`),
mirroring `docs/design.md` section 7:

### `raw_reanalysis`

FASTQ in, full pipeline runs:

```
FASTQC ─┐
        ├─> SAMPLE_QC (MultiQC)
TRIM_FASTP ─> SALMON_QUANT ─┘
                 │
        (SALMON_INDEX only if --salmon_index not supplied)
                 │
                 v
   DIFFERENTIAL_EXPRESSION_DESEQ2 (R/DESeq2, tximport)
                 │
                 v
          STANDARDIZE_OUTPUT
                 │
        ┌────────┴────────┐
        v                 v
  EMIT_MANIFEST      RENDER_REPORT
```

Required params: `--reads` (a `Channel.fromFilePairs` glob, e.g.
`'data/raw/GIGAS_HEAT01/*_R{1,2}.fastq.gz'`), `--gtf`, `--tx2gene`,
`--sample_sheet`, and either `--salmon_index` or `--transcript_fasta` (to
build one). See `nextflow.config` for the full parameter list and defaults
(all `null` — nothing is silently assumed).

### `processed_results_harmonization`

An existing DE-like TSV (e.g. already-published DESeq2/edgeR/limma output)
is validated and reshaped directly — no FASTQ processing is attempted:

```
(processed TSV) ─> STANDARDIZE_OUTPUT ─┬─> EMIT_MANIFEST (warns: raw QC unverified)
                                        └─> RENDER_REPORT
```

Required params: `params.rnaseq.processed_results`, `params.rnaseq.study_id`,
`params.rnaseq.comparison_id` — exactly the structure already defined in
`../../config/demo.config`, so this workflow is directly compatible with:

```bash
cd workflows/rnaseq
nextflow run main.nf -profile docker -config ../../config/demo.config
# or, without demo.config, the built-in convenience profile:
nextflow run main.nf -profile docker,demo

# On the lab workstation, add the `workstation` profile to keep the Nextflow
# work directory on the external SSD (see docs/first_raw_reanalysis.md,
# "Storage layout"):
nextflow run main.nf -profile local,workstation -c ../../config/CALLA2026_OSHV.config \
  --comparison_id midori_oshv1_france_vs_control
```

In this mode, `EMIT_MANIFEST` always records a warning in the manifest's
`warnings` array: raw QC (read counts, mapping rate, adapter content) could
not be independently verified, because no raw data passed through this
workflow run. This is the concrete implementation of the "quality_flags"
concept from `docs/design.md` section 7 at the workflow-manifest level (the
harmonize step in `src/harmonize/rnaseq.py` separately derives its own
evidence-record-level `quality_flags` from study/comparison metadata).

## Inputs and outputs

| | raw_reanalysis | processed_results_harmonization |
|---|---|---|
| Input | paired FASTQ + genome/transcriptome + GTF + sample sheet | one DE-like TSV |
| QC | FastQC, fastp, MultiQC (real) | not run; flagged in manifest |
| Quantification | Salmon (pseudoalignment + `--geneMap`) | n/a |
| DE model | DESeq2 via tximport | n/a (table assumed already DE-modeled upstream) |
| Standardized output | `STANDARDIZE_OUTPUT` validates/reshapes DESeq2 output | `STANDARDIZE_OUTPUT` validates/reshapes the provided table |
| Sample exclusion | `--exclude_samples`, recorded in the manifest | n/a |
| Viral read QC | optional, `--viral_reference` | n/a |

Both modes converge on the same standardized TSV schema, written to
`results/rnaseq/standardized/<study_id>_<comparison_id>_dge_standardized.tsv`:

```
gene_id  baseMean  log2FoldChange  lfcSE  stat  pvalue  padj
```

This is byte-for-byte the same column set as
`data/demo/rnaseq/GIGAS_HEAT01_acute_heat_vs_control_dge_demo.tsv` and
exactly what `src/harmonize/rnaseq.py::harmonize_rnaseq()` expects to read.

Other outputs (all under `params.outdir`, default `results/rnaseq/`):

- `manifests/<study_id>_<comparison_id>_manifest.json` — provenance manifest
  (input checksums, parameters, software versions, qc_metrics placeholder,
  warnings). Every process additionally emits a `versions.yml` in the
  nf-core convention (`process name -> tool -> version`).
- `reports/<study_id>_<comparison_id>_rnaseq_report.html` — Quarto-rendered
  HTML report (QC summary table, DE table preview, volcano plot) built from
  `assets/report_template.qmd`.
- (`raw_reanalysis` only) `fastqc/`, `trimmed/`, `salmon/`, `multiqc/`,
  `deseq2/` — intermediate per-step outputs.

## Sample exclusion and viral read QC (`raw_reanalysis` only)

**`--exclude_samples`** takes a comma-separated list of `sample_id`s to drop
before DESeq2. A name that is not in the sample sheet stops the run rather than
excluding nothing. The run's manifest records the list under
`parameters.exclude_samples` and adds a `samples_excluded_by_curation` warning;
the reason and supporting evidence belong in the study YAML
(`comparisons[].excluded_samples`), which `aree validate-study` checks against
the sample sheet. Exclude on an independent measurement — a control carrying
the pathogen, a sample swap — never because a sample weakens the result.

**`--viral_reference`** (a viral genome FASTA) turns on per-sample viral read
QC for pathogen-challenge studies. `VIRAL_QUANT` maps every trimmed read pair
of every library to the viral genome with Salmon; `VIRAL_SUMMARY` joins the
counts to the sample sheet and writes
`viral_load/<study_id>_<comparison_id>_viral_load.{tsv,json}`. Excluded samples
are still counted, and flagged `excluded_from_de`. The manifest carries the
JSON under `viral_load`, the report shows the table, and two warnings are
raised against `--viral_detection_per_million` (default 100):

- `viral_reads_in_control` — a control sample at or above the threshold;
- `no_viral_reads_in_challenged_sample` — a challenged sample below it.

The count is QC, not a titre. It scales with library depth and composition,
and a reference that differs from the challenge strain lowers counts in
infected animals; it cannot put reads into an uninfected one. For OsHV-1,
`config/CALLA2026_OSHV.config` points at RefSeq `NC_005881.2`
(`data/reference/pathogens/OsHV-1/`, with a provenance sidecar).

## Feeding output into the AREE CLI

The standardized TSV produced by either mode is the direct input to the
harmonization CLI step described in `CLAUDE.md`:

```bash
aree harmonize --study GIGAS_HEAT01 \
    --input results/rnaseq/standardized/GIGAS_HEAT01_acute_heat_vs_control_dge_standardized.tsv
```

which calls `src/harmonize/rnaseq.py::harmonize_rnaseq()` to convert the
table into evidence records (`schemas/evidence.schema.json`) joined against
the study's registry entry (`registry/studies/GIGAS_HEAT01.yaml`) for
tissue/life-stage/stressor/phenotype/sample-size context. The workflow
manifest JSON is a separate, complementary artifact — it documents how the
TSV was produced (or, in harmonization mode, that it was *not* independently
regenerated from raw data), not the evidence records themselves.

## What has and has not been run

- **Run**: the full `raw_reanalysis` path at full depth with `-profile local`
  on four real comparisons (three from PRJNA1329250, one from PRJNA593309), and
  `processed_results_harmonization` against the demo table. The 2026-09-30
  `midori_oshv1_australia_vs_control` re-run, with `--exclude_samples` and
  `--viral_reference`, is the first full-depth run to finish `OK`; earlier
  ones exited `ERR` at `RENDER_REPORT` on an exFAT work directory.
- **Not run**: any container from `../../containers/README.md`; the
  `-profile docker` or `-profile apptainer` paths. The non-stub run on the tiny
  CI fixtures fails in `SALMON_QUANT` under Salmon 2.6 ("failed to fill whole
  buffer" reading the trimmed fixture FASTQ); CI only stub-runs them.
- **Real, inspectable**: every `process` block declares a real upstream
  container tag, a real `label` (`process_low|process_medium|process_high`
  per `../../config/base.config`), typed `input:`/`output:` channels, and a
  `script:` block containing the actual command line (or, for the R and
  Python steps, an embedded real script) a bioinformatician would run. The
  DAG in `main.nf` uses real `Channel.fromFilePairs`/`Channel.fromPath`
  factories and real `if (params.mode == ...)` branching — it is not a
  linear echo-only stub chain.
- **CI coverage**: tiny paired FASTQ, reference, annotation, tx2gene, and sample
  fixtures live under `../../tests/fixtures/rnaseq_raw`. Every RNA-seq process
  has an explicit `stub:` output contract, so current Nextflow can exercise the
  complete raw DAG without a bioinformatics toolchain. This verifies channels,
  parameter validation, and artifact handoffs only; it is not a container test
  and produces no biologically meaningful result. Image tags in
  `../../containers/README.md` remain unverified.

## Profiles

Defined in `nextflow.config`:

- `docker` — runs all containers via Docker (`docker.enabled = true`).
- `apptainer` — runs all containers via Apptainer/Singularity, pulling the
  same `docker://` image references (for HPC use; untested here).
- `local` — disables container engines entirely (assumes tools are on
  `$PATH`; not recommended, provided for completeness).
- `demo` — convenience profile that points
  `params.rnaseq.{processed_results,study_id,comparison_id}` at the
  repository's synthetic heat-stress demo TSV, equivalent to running with
  `-config ../../config/demo.config`.

Combine an execution profile with `demo`, e.g. `-profile docker,demo`.
