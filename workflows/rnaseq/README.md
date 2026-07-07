# AREE :: RNA-seq workflow

Implements CLAUDE.md Layer 2, section A (RNA-seq): FASTQ QC, adapter
trimming, pseudoalignment/quantification, sample QC, differential
expression, effect-size calculation, and a standardized, enrichment-ready
result table plus a human-readable report.

**Status: structurally complete Nextflow DSL2 scaffold. Nothing in this
directory has been executed against real data in this build.** There was no
compute budget and no real (or synthetic) FASTQ/BAM data available in this
environment. Every process declares a real container image and a real tool
command line, and the DESeq2 step embeds a real, syntactically valid R
script — but no container has been pulled, no `nextflow run` has been
invoked, and no output file in this repository was produced by actually
running this pipeline. Treat everything here as reviewed-but-unexecuted
scaffold code, consistent with `docs/design.md` section 9 ("Explicit
assumptions").

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

## What has and has not been run in this build

- **Not run**: `nextflow run main.nf` has not been invoked. No container in
  `../../containers/README.md`'s table has been pulled or executed. No
  FASTQ, BAM, Salmon index, or DESeq2 RData file in this repository was
  produced by this pipeline.
- **Real, inspectable**: every `process` block declares a real upstream
  container tag, a real `label` (`process_low|process_medium|process_high`
  per `../../config/base.config`), typed `input:`/`output:` channels, and a
  `script:` block containing the actual command line (or, for the R and
  Python steps, an embedded real script) a bioinformatician would run. The
  DAG in `main.nf` uses real `Channel.fromFilePairs`/`Channel.fromPath`
  factories and real `if (params.mode == ...)` branching — it is not a
  linear echo-only stub chain.
- **Known gaps before this could actually run**: no synthetic FASTQ/BAM
  fixtures exist anywhere in this repository (see `../../config/demo.config`
  header), so `raw_reanalysis` mode cannot currently be smoke-tested even in
  CI; only `processed_results_harmonization` mode is demo-able today. Image
  tags in `../../containers/README.md` have not been re-verified as current.
  The `DIFFERENTIAL_EXPRESSION_DESEQ2` process includes a `stub:` block
  (usable with `nextflow run -stub-run`) so the DAG can at least be
  structurally exercised without a real Bioconductor container.

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
