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

**Single-end libraries.** Add `--single_end` and point `--reads` at one file
per library, e.g. `'data/raw/DUPERRET2025_STARVE_RNASEQ/ERR*.fastq.gz'`. The
sample id is the file name without `.fastq`/`.fq`(`.gz`). For ENA downloads
that is the run accession, which the DESeq2 step matches to the sample
sheet's `run_accession` column. The single-end path changes three steps:

- fastp trims a **given** 3′ adapter, `--se_adapter` (default: Illumina
  TruSeq/NEBNext read 1, `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`). It does not
  auto-detect, because fastp's single-end detection guesses from
  overrepresented sequence. On real DUPERRET2025_STARVE_RNASEQ reads it chose
  an abundant oyster transcript and trimmed it from about 1% of reads, with
  a separate guess per library. Pass `--se_adapter auto` to restore detection
  deliberately.
- Salmon quantifies with `-r` and without `--gcBias`. Fragment GC-bias
  correction needs the fragment-length distribution that single-end reads do
  not observe, and Salmon marks it beta for single-end data.
- Viral read QC counts reads; its `*_pairs` columns hold reads.

`single_end` and `se_adapter` are recorded in the manifest's `parameters`.

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
| Covariates | `--covariates` (e.g. `family`), fitted design recorded in the manifest | n/a |
| Gene exclusion | `--exclude_genes` (e.g. rRNA in total-RNA libraries), counts recorded in the manifest | n/a |
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

## 3' Tag-seq libraries (`raw_reanalysis` only)

`--tag_seq` is for 3' tag libraries such as Lexogen QuantSeq, which give
about one read per transcript whatever its length. It is independent of
`--single_end` (QuantSeq also comes paired-end) and changes three steps:

| Step | Standard RNA-seq | `--tag_seq` |
|---|---|---|
| fastp | adapter trimming | also `--trim_poly_x` (poly-A tails) |
| Salmon | `--gcBias` (paired-end) | `--noLengthCorrection`, no `--gcBias` |
| DESeq2 | `DESeqDataSetFromTximport` (per-sample length offsets) | `DESeqDataSetFromMatrix` on the counts, no length offset; pseudo-bulk sums raw counts |

Applying length correction to 3' tags makes long genes look less expressed
than they are. A quick way to tell a QuantSeq FWD library from full-length
RNA-seq when the deposit does not say: the fraction of reads containing 15+
consecutive A's is several percent for QuantSeq and about 0.1% for
full-length libraries. `tag_seq` is recorded in the workflow manifest's
`parameters`, and the length model DESeq2 actually used in `deseq2_design`
(`length_model`).

The length offsets do more than shift expression levels. They can create fold
changes. Measured on simulated QuantSeq-like reads from the *M.
galloprovincialis* RefSeq transcriptome (six single-end libraries, 3 vs 3,
2026-10-12): both modes recovered 59 of 60 genes with a simulated 4-fold
change. But for 40 genes whose expressed isoform switched to one about 2.7
times longer, with no change in 3' reads, the standard path called 26
differentially expressed (median |log2FC| 0.84, up to 3.0). The per-sample
length offset moves with the isoform mix even when the 3' read count does not.
With `--tag_seq`, none were called (median |log2FC| 0.09).

## Covariates (`raw_reanalysis` only)

**`--covariates`** takes a comma-separated list of sample-sheet columns to fit
as factors ahead of the contrast, e.g. `--covariates family` fits
`~ family + condition`. The contrast is still condition. Use it for a known
source of variation that is balanced across both arms, such as family,
sequencing batch or lane. The run stops, naming the problem, when:

- a column is missing from the sample sheet, or empty for a sample in the contrast;
- a covariate is `condition`, an identifier, or the `--replicate_unit` column;
- a covariate is confounded with condition (the model is not full rank),
  e.g. one family per arm;
- the model leaves no residual degrees of freedom;
- with `--replicate_unit`, a covariate varies within a unit (summed units
  cannot carry it).

A covariate with only one level in the contrast is dropped, not fitted. This
happens when a run is subset to one family. The drop is reported as a
`covariates_dropped_constant` warning.

DESeq2 writes `deseq2/<study_id>_<comparison_id>_deseq2_design.json` (the
formula fitted, covariates requested/used/dropped, libraries or units, and
residual df). The manifest carries it under `deseq2_design` and the request
under `parameters.covariates`.

**A covariate does not fix clustering.** If several libraries share the unit
that received the treatment, they are not independent replicates, whatever
covariates are fitted. DESeq2 cannot take that unit as a term when it is
nested in condition. An example is DUPERRET2025_STARVE_RNASEQ, where both
families share each tank. There, `~ family + condition` would count 6
libraries per diet when 3 tanks were treated. The defensible fits are
`--replicate_unit tank` (sum each tank's libraries, 3 vs 3) or one run per
family, as the source paper did.

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

## Gene exclusion for total-RNA libraries (`raw_reanalysis` only)

**`--exclude_genes`** takes a file of `gene_id`s, one per line, with `#`
starting a comment. Those genes are removed from the tximport matrices before
anything else is estimated: pseudo-bulk sums, size factors, dispersions and
the contrast. The use it was written for is rRNA in libraries that were
neither rRNA-depleted nor poly-A selected. In `ZHONG2025_PONV1_MORTALITY`
(field metatranscriptomes), rRNA took 65–97% of gene counts, and its share
differed between farm and wild animals.

What it does and does not change was measured on that study. Three farm and
three wild mantle libraries were run at 1 M read pairs each, with and without
the rRNA list:

- **Size factors are unchanged** (every ratio 1.00), and so are fold
  changes for the other genes (r = 1.00). DESeq2's median-of-ratios takes
  a median over thousands of genes, and the six rRNA genes Salmon reports
  (it collapses the near-identical rRNA copies) cannot move it. Without the
  option, normalization is not biased.
- **rRNA genes stop appearing as results.** Four rRNA genes were in the
  uncorrected result table (499 vs 496 genes at padj < 0.05). They would
  otherwise reach the evidence table as "differentially expressed genes".
- **The rRNA share per library is recorded** (here 73–97%). It is the
  number to read when an arm's libraries look different.

Use it for any total-RNA study. It is a cleanliness and QC step, not a
correction to normalization.

Build the list from the annotation the run uses, so the ids match `tx2gene`:

```bash
python scripts/list_genes_by_biotype.py \
    --gtf data/reference/GCF_963853765.1/genomic.gtf --biotype rRNA \
    --tx2gene data/reference/GCF_963853765.1/tx2gene.tsv \
    --out data/reference/gene_lists/GCF_963853765.1-RS_2024_06_rRNA.txt
```

The list for the current reference is committed:
`data/reference/gene_lists/GCF_963853765.1-RS_2024_06_rRNA.txt`, 68 genes
with NCBI `gene_biotype "rRNA"`. Only those present in the quantification are
dropped: Salmon's index collapses identical sequences, so many rRNA copies
have no row of their own. The manifest's `n_excluded` counts the ones that
did. A list that matches no quantified gene stops
the run, so a list built for another annotation cannot silently do nothing.
The manifest records the path under `parameters.exclude_genes`. Under
`deseq2_design.exclude_genes` it records the file's name and MD5, how many
genes were listed and dropped, and the min, median and max share of each
library's counts they held. That share is QC worth reading: a large
difference between arms means the arms' libraries differ in composition, not
just in depth.

Excluded genes get no result rows, so they never reach the evidence table.
Use this for technical contaminants of the count matrix (rRNA,
mitochondrial rRNA, spike-ins), never to remove genes because of their
results.

## Several comparisons of one study

Run each registered comparison as its own `nextflow run`, one after another,
from the same launch directory, with the same work directory, and with
`-resume`. Give each run its own `--outdir`:

```bash
for c in starved_vs_fed_0hpc fed_oshv1_24hpc_vs_0hpc; do
  nextflow -c config/DUPERRET2025_STARVE_RNASEQ.config run workflows/rnaseq \
    -profile local,workstation -resume --comparison_id "$c" ... \
    --outdir "results/rnaseq/DUPERRET2025_STARVE_RNASEQ/$c"
done
```

The first run does all the per-library work. Later runs take FASTQC, fastp,
`SALMON_INDEX`, `SALMON_QUANT`, `VIRAL_INDEX`, `VIRAL_QUANT` and MultiQC from
the cache. They rerun only DESeq2, viral summary, standardization, manifest
and report. The per-library processes read study-level params only (`reads`,
`single_end`, `se_adapter`, `tag_seq`, the reference files), so a new
`--comparison_id`, `--control_level`/`--treatment_level`, `--replicate_unit`,
`--exclude_samples` or `--outdir` leaves their task hashes unchanged. A cached
task's outputs are published to the new `--outdir` as well.
`tests/test_rnaseq_resume_across_comparisons.py` checks this with two stub
runs, paired and single end, in CI's Nextflow job.

**Launch Nextflow with the `nextflow` launcher script**, e.g. from
`conda activate aree` (`environment.yml`) or a current launcher from
nextflow.io. Do not run `java -jar ~/.nextflow/framework/<ver>/nextflow-<ver>-one.jar`.
The launcher passes `--add-opens` flags to the JVM. Without them, Nextflow
cannot serialize task cache entries. The writes fail silently, every run
reports `cached=0`, and every library is recomputed. `main.nf` checks for
this at startup. It warns on a plain run and stops a `-resume` run.

That happened with the three DUPERRET2025_STARVE_RNASEQ comparisons on
2026-10-12. `-dump-hashes` showed the same FASTQC hash for a library in every
comparison. The cache database held no entries, so Nextflow treated each
existing task folder as foreign and moved to a new hash. That is why one
library showed different task folders in different runs. Those folders are not
in any cache. The first correctly launched run recomputes them once, and later
comparisons reuse its tasks. Remove the old work directory only after that run
has finished.

## Feeding output into the AREE CLI

The standardized TSV produced by either mode is the direct input to the
harmonization CLI step described in `CLAUDE.md`:

```bash
aree harmonize --study GIGAS_HEAT01 \
    --input results/rnaseq/standardized/GIGAS_HEAT01_acute_heat_vs_control_dge_standardized.tsv
```

which calls `src/harmonize/rnaseq.py::harmonize_rnaseq()` to convert the
table into evidence records (`schemas/evidence.schema.json`) joined against
the study's registry entry (`registry/demo_studies/GIGAS_HEAT01.yaml`) for
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
- **Single-end**: executed natively (`-profile local`, Salmon 2.3.4, fastp
  1.3.6) on 2026-10-11 against 40,000 real reads from each of four
  DUPERRET2025_STARVE_RNASEQ runs and a 300-transcript RefSeq subset. All 24
  tasks completed, through DESeq2, the manifest and the report. That checks
  the command lines, not the biology; no full-depth single-end run has been
  made.
- **3' Tag-seq**: executed natively on 2026-10-12, with and without
  `--tag_seq`, against six simulated QuantSeq-like single-end libraries
  (82,000 to 101,000 reads) and the full *M. galloprovincialis* RefSeq
  transcriptome; see the Tag-seq section for the result. Salmon received
  `--noLengthCorrection` and fastp `--trim_poly_x` in every library. No real
  Tag-seq study has been run.
- **CI coverage**: tiny paired FASTQ, reference, annotation, tx2gene, and sample
  fixtures live under `../../tests/fixtures/rnaseq_raw`, and single-end ones
  (ENA-style file names, `run_accession` in the sample sheet) under
  `../../tests/fixtures/rnaseq_raw_se`. Every RNA-seq process
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
