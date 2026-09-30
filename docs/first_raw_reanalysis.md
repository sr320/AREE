# The first raw-reanalysis study

`CALLA2026_OSHV` (BioProject `PRJNA1329250`) is the second real study
registered in AREE and the first in `raw_reanalysis` mode. This page records
what was curated, the judgement calls made, and — plainly — what has not been
run.

**Status: three of six comparisons run at full depth and harmonized; one pools.**
`qc_status` and `analysis_status` are both `in_progress`.

| Comparison | State | Evidence records | In meta-analysis |
|---|---|---|---|
| `miyagi_oshv1_usa_vs_control` | full depth, harmonized | 30,625 | **yes** — the prespecified `meta_analysis_primary` comparison |
| `midori_oshv1_australia_vs_control` | full depth, harmonized, DO_CT_4 excluded (10 samples) | 30,482 | no — within-study, not primary |
| `midori_oshv1_france_vs_control` | full depth, harmonized, DO_CT_4 excluded (10 samples); only 2 of 5 challenged animals infected | 30,541 | no — within-study, not primary |
| `midori_oshv1_usa_vs_control` | not run; DO_CT_4 exclusion already registered | 0 | no |
| `miyagi_oshv1_australia_vs_control` | not run | 0 | no |
| `miyagi_oshv1_france_vs_control` | not run | 0 | no |

The Miyagi USA comparison pools with `DELISLE2020_OSHV_TEMP` (`PRJNA593309`)
across 23,094 shared genes — see [Full-depth runs](#full-depth-runs). Per-sample
viral read QC found an infected control and three uninfected challenged
animals — see [Viral read QC](#viral-read-qc-and-the-do_ct_4-exclusion). The
pipeline was first proven on a subsampled slice; that history is kept in
[What the pilot found](#what-the-pilot-found).

## Why this study, and why raw data

`HESSER2024_VCOR` exposed the limit of harmonizing published tables: it
resolves 87.2% of its identifiers, harmonizes cleanly, and contributes
**nothing** to any pooled estimate, because its supplementary table reports
only `log2FoldChange` and `padj`. No standard error, no inverse-variance
weight, no pooling.

Reanalyzing raw reads removes that dependency. AREE's own DESeq2 run emits
`lfcSE` and an unadjusted p-value by construction, whatever the authors chose
to publish. That is the reasoning behind the selection criterion in
[candidate_studies.md](candidate_studies.md), and this study was picked to test
it: 42 paired-end libraries, all groups replicated, and the same stressor class
as the study already registered, so the two can eventually pool.

## The study

Calla B, Thompson NF, Burge CA (2026). *Population-specific transcriptomics of
Pacific oyster after exposure to a highly pathogenic, globally distributed
virus.* Fish & Shellfish Immunology 171:111154.
[10.1016/j.fsi.2026.111154](https://doi.org/10.1016/j.fsi.2026.111154)

USDA-ARS Pacific Shellfish Research Unit, Newport, Oregon. Spat from two
hatchery populations challenged with three OsHV-1 microvariant isolates.

| | Control | Australia | France | USA |
|---|---|---|---|---|
| **Midori** | 6 | 5 | 5 | 5 |
| **Miyagi** | 6 | 5 | 5 | 5 |

42 runs · Illumina NovaSeq 6000 · paired-end · whole soft tissue · 226 GB.

The design was **not transcribed from the paper**, which is paywalled. It was
read programmatically from the deposited BioSample attributes:

```bash
aree fetch-samplesheet --bioproject PRJNA1329250 --study CALLA2026_OSHV \
  --condition-attribute breed --condition-attribute "Viral strain" \
  --attribute dev_stage --attribute tissue
```

That writes `samplesheet.tsv`, `fastq_manifest.tsv` (with ENA's own MD5 per
file), and `ena_provenance.json` recording the queries and checksums. Reading
the design from the archive rather than the methods section is not a
convenience — it is how an unreplicated design gets caught before anyone
downloads 226 GB. `aree validate-study` cross-checks that the study's declared
BioProject matches the one the sample sheet was generated from.

## Curation judgements

**This is not registered as resilience evidence.** The BioProject is titled
*"Evaluating Pacific oyster lineages for tolerance to Ostreid herpesvirus"*,
and it would have been easy to register the comparisons as
`disease_resistance`. Neither the deposited metadata nor the publication
abstract reports a survival, mortality, or viral-load measurement for these
animals; the paper's own framing is that its results matter for *future* work
developing tolerant oysters. The comparisons are therefore registered as
`disease_susceptibility` / `disease` and flagged
`ambiguous_phenotype_definition`. A population difference in transcriptional
response is not a tolerance phenotype — see
[resilience_vs_exposure.md](resilience_vs_exposure.md).

**Unknown exposure parameters are null, not guessed.** Dose, route, duration,
and sampling timepoint are absent from the deposited metadata and the paper is
not open access. `exposure_intensity`, `exposure_duration`, and
`exposure_timing` are empty, and a test asserts they stay that way. Effect
sizes from this study are therefore **not dose-comparable** to another
pathogen-challenge study until someone reads the methods and fills them in.

**The three viral isolates are kept as separate comparisons.** The publication
reports no marked difference between the USA and Australian variant responses.
That is the authors' conclusion about their data, not a licence for AREE to
merge the arms; pooling them is a decision for meta-analysis, with heterogeneity
reported, not a curation shortcut.

**Population spelling.** The deposited BioSample records spell the second
population `Myagi`; the publication spells it `Miyagi`. The deposited spelling
is preserved verbatim in the sample sheet, and the discrepancy is recorded in
`strain_or_population`.

**Reference choice.** AREE reanalyzes against `GCF_963853765.1` /
`xbMagGiga1.1`, annotation release `RS_2024_06` — the same annotation the
identifier crosswalk is built from. This study therefore carries **no crossing**
between its assembly and AREE's, unlike `HESSER2024_VCOR`
([handling_genome_versions.md](handling_genome_versions.md)). The paper's
keywords mention a transcriptome assembly, so the authors' own analysis may
have used a different reference: AREE's reanalysis is an independent
reanalysis, not a replication, and is not expected to reproduce their figures.

## Running it

The commands below are the ones actually used for the pilot, adapted to full
depth. The pilot ran natively rather than in containers; a `-profile docker`
run has **not** been tried.

```bash
# Toolchain (Apple Silicon; no Docker required).
brew install nextflow fastp fastqc salmon
Rscript -e 'BiocManager::install(c("DESeq2","tximport"))'
pip install multiqc

# 1. Stage the reads (226 GB at full depth), then verify ENA's checksums.
awk -F'\t' 'NR>1 {print $4"  "$1"_"$2".fastq.gz"}' \
  data/studies/CALLA2026_OSHV/fastq_manifest.tsv > md5sums.txt
md5sum -c md5sums.txt

# 2. Reference transcriptome + annotation, and a transcript->gene map.
BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/853/765/GCF_963853765.1_xbMagGiga1.1
curl -O "$BASE/GCF_963853765.1_xbMagGiga1.1_rna.fna.gz"
curl -O "$BASE/GCF_963853765.1_xbMagGiga1.1_genomic.gtf.gz"
awk -F'\t' '$3=="transcript"{
  if (match($9, /transcript_id "[^"]+"/)) tx=substr($9,RSTART+15,RLENGTH-16)
  if (match($9, /GeneID:[0-9]+/))        g=substr($9,RSTART+7,RLENGTH-7)
  if (tx && g) print tx"\t"g
}' genomic.gtf | sort -u > tx2gene.tsv

# 3. Run one comparison. The `workstation` profile puts the work directory on
#    the T5 SSD (see "Storage layout" below); never put it on exFAT.
nextflow run workflows/rnaseq -profile local,workstation \
  -c config/CALLA2026_OSHV.config \
  --control_level Midori_Control \
  --treatment_level Midori_France \
  --comparison_id midori_oshv1_france_vs_control

# 4. Harmonize the standardized output. --comparison is required: the workflow
#    names its output for the stage that produced it, and a raw-mode study has
#    no results_file to match a filename against.
export AREE_CROSSWALK=data/reference/crosswalk/mgigas_gene_id_crosswalk.tsv
aree harmonize --study CALLA2026_OSHV \
  --comparison midori_oshv1_france_vs_control \
  --input results/rnaseq/standardized/CALLA2026_OSHV_midori_oshv1_france_vs_control_dge_standardized.tsv
```

Then commit the standardized table, workflow manifest, and report under
`data/studies/CALLA2026_OSHV/`, set that comparison's `results_file` in the
study YAML, and keep `analysis_status` current — `aree harmonize` warns if it
still reads `not_started`. Do **not** set `meta_analysis_primary` on a new
comparison: this study already has its one prespecified comparison, and a
second flag brings back the within-study error described below.

Harmonize one comparison at a time with `--comparison`, as above. Plain
`aree harmonize --study CALLA2026_OSHV` stops at the first comparison whose
`results_file` is still null, and four are.

## What the pilot found

Before committing to a 226 GB download, one contrast
(`midori_oshv1_france_vs_control`, 11 libraries) was run at reduced depth: the
first 1,000,000 read pairs per library, ~1.5 GB total, streamed from ENA and
truncated. The purpose was to make the workflow run, not to produce biology.

It was the right order of operations. **The pipeline could not have worked**,
and a full-depth run would have died in the first seconds.

### The run

38 processes, all green: FASTQC ×22 → fastp ×11 → Salmon index → Salmon quant
×11 → MultiQC → DESeq2 → standardize → manifest → report. Roughly 7 minutes on
an M4 with 10 cores, run natively via `-profile local` with Homebrew-installed
tools; no container was pulled.

Two results are worth recording:

* **22,301 genes quantified, every one with a real `lfcSE`**, and 22,276 with an
  unadjusted p-value. This is precisely what `HESSER2024_VCOR` cannot supply and
  the entire reason this study was selected for raw reanalysis.
* **Identifier resolution is 100.0% exact** (4,000-identifier dry run against
  the real crosswalk), versus 87.2% for `HESSER2024_VCOR`. This confirms the
  prediction recorded in the study YAML: quantifying against the same annotation
  the crosswalk is built from removes the cross-assembly step entirely. All
  33,068 GeneIDs in the reference GTF are present in the crosswalk, none
  unmatched.

**The pilot output is not valid biology and has not been harmonized.** Taking
the first N reads of a FASTQ is not a random subsample — early reads come from
one region of the flowcell — and 1M pairs badly under-powers a differential
expression test. The numbers above describe pipeline behaviour, nothing else.

### Sixteen defects

| # | Where | Defect |
|---|---|---|
| 1 | rnaseq, proteomics | bare script-level statements (`def VALID_MODES`, `workflow_version`) — **would not compile** on Nextflow 26.04 |
| 2 | rnaseq, proteomics, methylation | `workflow.onComplete` at script level — **would not compile** |
| 3 | same | `params` unresolvable inside the relocated closure — NPE after every run |
| 4 | same | `workflow` unresolvable in the same closure — NPE after every run |
| 5 | metabolomics | `Channel.empty().collect()` emits nothing, so `EMIT_MANIFEST` and `RENDER_REPORT` were **silently skipped**, exit 0 |
| 6 | proteomics | `RENDER_REPORT` copied the template onto itself; `cp` rejects identical paths |
| 7 | DESeq2 module | condition levels hardcoded to `c("control","treatment")` — breaks **every** real study |
| 8 | DESeq2 module | `quant_subdir` mandatory though Salmon names dirs for the sample |
| 9 | DESeq2 module | `tx2gene` read with `header=TRUE` unconditionally, silently dropping the first transcript |
| 10 | `aree harmonize` | comparison resolved by filename only; workflow output matched nothing |
| 11 | `aree harmonize` | **`--input` was ignored** — it selected a comparison, then harmonized the registry's `results_file` instead of the file passed |
| 12 | `aree harmonize` | manifest writing crashed on any path outside the repo |
| 13 | SAMPLE_QC | declared `multiqc_data`; MultiQC writes `multiqc_report_data`, so the task failed as a missing output after exiting 0 |
| 14 | DESeq2 module | `ignoreTxVersion` strips versions from quant.sf but not from `tx2gene`, so any NCBI GTF stops matching |
| 15 | DESeq2 module | a bare `$` before a quote is mangled by Groovy interpolation, corrupting the R regex |
| 16 | rnaseq banner | printed `comparison_id: null` for a value the run was using |

Defect 11 deserves emphasis. `aree harmonize --study X --input results.tsv` is
one of the six commands specified in the project brief, and it never harmonized
the file it was given. The existing tests missed it because they always passed
the path that was already the declared `results_file`, so reading the wrong file
produced identical output. It also made raw reanalysis impossible by
construction: a raw-mode study has `results_file: null`, so there was nothing to
fall back to.

Defect 5 is the most dangerous kind: a run that exits 0, reports
`completed=1`, and silently produces no manifest and no report. Every other
defect on this list announces itself.

### Storage layout

The internal disk is too small for raw reads, so the repository reaches
external storage through two symlinks and one Nextflow profile. Both drives are
HFS+ or APFS, which support the symlinks and permissions Nextflow relies on.

| Location | Drive | Purpose |
|---|---|---|
| `/Volumes/Samsung T5/aree/work` | Samsung T5, SSD, ~450 MB/s | Nextflow work directory (`-profile workstation`) |
| `/Volumes/blue-block/aree/raw/<STUDY_ID>/` | blue-block, HDD, 4 TB, ~140 MB/s | Raw FASTQs; `data/raw` is a symlink here |
| `/Volumes/blue-block/aree/reference/<ASSEMBLY>/` | blue-block | Transcriptome, GTF, tx2gene, Salmon index; `data/reference/<ASSEMBLY>` is a symlink here |
| `results/` | internal APFS | Published outputs, unchanged, so `aree harmonize` relative paths keep working |

The symlinks are gitignored (`data/raw/`, `data/reference/GCF_*/`). Recreate
them on a new machine with:

```bash
ln -sfn /Volumes/blue-block/aree/raw data/raw
ln -sfn /Volumes/blue-block/aree/reference/GCF_963853765.1 data/reference/GCF_963853765.1
```

The `workstation` profile lives in `config/base.config`, so every assay
workflow accepts it. Combine it with an executor profile
(`-profile local,workstation`); it sets only `workDir`.

### Environment notes

* **Quarto fails on exFAT.** `RENDER_REPORT` dies in Quarto's cleanup
  (`safeRemoveSync`) when the Nextflow work directory is on an exFAT volume;
  the identical run succeeds on APFS. Publishing Salmon output directories to
  exFAT also fails. Keep `-work-dir` and `--outdir` on APFS even when read data
  lives on an external drive.
* **Do not begin a read glob with a character class.** `[SD]*_{1,2}.fastq.gz`,
  used to skip the `._` AppleDouble sidecars macOS creates on exFAT, silently
  breaks `fromFilePairs` key extraction and yields empty sample ids. Anchor on
  the run-accession prefix (`SRR*`) instead.
* No container was used. `-profile local` with `brew install nextflow fastp
  fastqc salmon` plus `BiocManager::install(c("DESeq2","tximport"))` and
  `pip install multiqc` is sufficient on Apple Silicon, and avoids x86
  emulation. The container tags in `containers/README.md` remain unverified.

## Full-depth runs

Three comparisons have been run on all reads and committed, 11 libraries each
(6 control, 5 challenged), with the same `-profile local` toolchain as the
pilot, against `GCF_963853765.1` / `RS_2024_06`. Each committed workflow
manifest records the input checksum and parameters. Its `software_versions`
block lists the versions the workflow declares for its containers (Salmon
1.10.3, DESeq2 1.42.0, tximport 1.30.0), not the Homebrew builds that actually
ran — Salmon 2.6.0 (per each quant directory's `cmd_info.json`) and, at the
2026-09-30 re-run, DESeq2 1.50.2 and tximport 1.38.2. Until the manifest
records observed versions, treat that block as intent rather than record.

The table below is the state before the DO_CT_4 exclusion. Re-run without
DO_CT_4 on 2026-09-30, Midori × Australia has 30,482 genes, 5,473 at
`padj < 0.05`, and a median |log2FC| of 0.38; Midori × France has 30,541
genes, 2,015 at unadjusted `p < 0.05`, 253 at `padj < 0.05`, and a median
|log2FC| of 0.21. See [Viral read QC](#viral-read-qc-and-the-do_ct_4-exclusion).

| | Miyagi × USA | Midori × Australia | Midori × France |
|---|---|---|---|
| Workflow start | 2026-08-31 | 2026-09-01 | 2026-09-04 |
| Genes in the standardized table | 30,625 | 30,561 | 30,628 |
| With `lfcSE` | 30,625 | 30,561 | 30,628 |
| With an unadjusted p-value | 30,522 | 30,376 | 30,527 |
| Unadjusted `p < 0.05` | 7,428 | 8,053 | 1,659 |
| DESeq2 `padj < 0.05` | 4,236 | 5,074 | **54** |
| Median \|log2FC\| | 0.31 | 0.35 | 0.19 |
| Identifier mapping | 100% `exact` | 100% `exact` | 100% `exact` |

**Midori × France shows almost no response.** With 1,659 genes at unadjusted
`p < 0.05` out of ~30,500 tested (about 1,500 expected by chance alone) and 54
after adjustment, this contrast is close to null. Viral read QC explains most
of it: only two of the five challenged animals carry the virus (next section).
Do not read the near-null result as Midori tolerance.

The pilot's prediction held at full depth: quantifying against the crosswalk's
own annotation leaves no identifier unresolved.

**Every full-depth run so far ended with Nextflow status `ERR`**, and in each
case the only failed task was `RENDER_REPORT`, the known Quarto cleanup failure
on an exFAT work directory (the work directory was on the external `Zarra`
volume, against the advice in [Environment notes](#environment-notes)). Every
task up to and including `STANDARDIZE_OUTPUT` and `EMIT_MANIFEST` exited 0.
Quarto writes the HTML before its cleanup step crashes, so a report does exist
in the failed task's directory. For Miyagi USA, the committed report was instead
rendered by hand from an APFS staging directory.

Each run's `--outdir` was under `/private/tmp`, which does not survive a
reboot. The Midori France outputs were lost from there and recovered on
2026-09-24 from the Nextflow work directory on `Zarra`: the standardized table
from `STANDARDIZE_OUTPUT` (`work/3a/6814ae…`), the manifest from `EMIT_MANIFEST`
(`work/86/108a74…`), and the report from the failed `RENDER_REPORT`
(`work/c4/ba4ef4…`). The table's SHA-256 matches the manifest. Point
`--outdir` at persistent storage for the remaining runs. Both manifests leave `qc_metrics` null (`n_samples`, read
depth, mapping rate): the workflow does not yet carry Salmon and MultiQC figures
into the manifest, so per-library QC lives only in the MultiQC output, which is
not committed.

### Viral read QC and the DO_CT_4 exclusion

The deposited metadata record no viral load, so AREE measures it: the RNA-seq
workflow's optional `--viral_reference` step maps every read pair of every
library to the OsHV-1 genome (RefSeq `NC_005881.2`,
`data/reference/pathogens/OsHV-1/`) and reports viral pairs per million
([workflows/rnaseq/README.md](../workflows/rnaseq/README.md)). It was first run
by hand on 2M-pair subsamples (2026-09-24), then on all reads by the workflow
(2026-09-30).

| Animals | OsHV-1 pairs per million |
|---|---|
| Miyagi controls (6), subsample | 0–1.5 |
| Miyagi × USA (5 of 5), subsample | 8,565–23,533 |
| Midori controls CT_1, 2, 3, 5, 6, all reads | 0–0.24 |
| **Midori control DO_CT_4, all reads** | **11,597** |
| Midori × Australia (5 of 5), all reads | 14,269–80,784 |
| Midori × France FR_2, FR_3, all reads | 11,096 and 57,978 |
| Midori × France FR_1, FR_4, FR_5, all reads | 0, 38.6, 0.3 |

Two findings follow, and they agree with the host transcriptome: the same
animals separate on the first principal component, and ADAR (`LOC105341503`),
an OsHV-1 response gene, sits at log2 CPM 7.6–9.5 in exactly the virus-positive
animals and 2.3–3.4 in the rest.

**DO_CT_4 is not a valid control.** An unchallenged animal carrying the virus at
infected-animal levels is either infected or mislabeled; the data cannot say
which. It is shared by every Midori comparison, so it is excluded from all
three (`excluded_samples` in the study YAML; `--exclude_samples DO_CT_4` at run
time). The exclusion rests on the viral measurement, not on its effect on the
result — but the effect is large. Without it, Midori × Australia gains 399
genes at `padj < 0.05`, and ADAR moves from log2FC 2.05 (padj 0.27, not
significant) to 5.99 (padj 1.6 × 10⁻⁷²): the infected control was masking the
strongest antiviral signal while leaving the genome-wide fold changes largely
unchanged (r = 0.97). Midori × France shows the same pattern at smaller
scale: 54 → 253 genes at `padj < 0.05`, and ADAR from log2FC 1.47 (padj 1.0)
to 5.42 (padj 0.003).

**Midori × France: 2 of 5 challenged animals infected.** All five are kept.
Dropping uninfected animals after seeing the data would turn "response to
challenge" into "response to infection" — a different comparison, which could
be registered separately if wanted. Whether FR_1, FR_4, and FR_5 resisted
infection or were never effectively exposed cannot be told from reads alone.

The viral count is QC, not a titre: it scales with depth, and the reference is
the original OsHV-1 genome rather than the microvariants used here, which
lowers counts in infected animals but cannot put reads into an uninfected one.
Per-animal viral load is nonetheless the closest thing to a measured phenotype
this study has.

### Why only one comparison pools

Two comparisons from one study cannot enter the same pool as if they were
independent studies. The Midori and Miyagi arms use different animals and
different controls, but they come from one lab and one challenge
experiment. Counting them separately would also let this study, on its own,
meet the "two independent studies" bar for a `high_priority_cross_study`
candidate. `aree meta-analyze` therefore refuses a group that holds two
comparisons from one study.

The way out is for the curator to prespecify one comparison per study with
`meta_analysis_primary: true`. Here that is `miyagi_oshv1_usa_vs_control`,
chosen on 2026-09-24 because it was the first run at full depth and already
formed the OsHV-1 pool. The Midori Australia effects were not compared against
it before the choice was made, but the Miyagi USA results had already been
examined, so the choice was not blind to them. The rationale is recorded in the
study's `curation_notes`; Midori France was recovered after the choice was
made. Records from both Midori comparisons stay in the evidence table and on
evidence cards; the meta-analysis counts them in
`n_excluded_non_primary`.

Note also that the completed comparisons never pair **both populations with
the same isolate** (Miyagi has only USA; Midori has Australia and France), so
they do not give a like-for-like Midori-versus-Miyagi contrast. That needs the two populations challenged with the same isolate — for
example `midori_oshv1_usa_vs_control` next to the Miyagi USA run.

### The pool

With the Miyagi USA comparison and `DELISLE2020_OSHV_TEMP`
(`oshv1_21c_96h_vs_21c_0h`), `aree meta-analyze --phenotype
disease_susceptibility --feature-type gene` produces 30,720 real-data features:

* **23,094 pooled across both studies (k = 2)**, of which 2,123 have a
  Benjamini–Hochberg `adjusted_p_value < 0.05` within the 30,720-test family;
* 4,387 of the two-study features have I² > 50%, so between-study heterogeneity
  is common and pooled effects should be read with it;
* the rest are single-study features, with nothing to replicate them.

Every one of these is `disease_associated` evidence. Neither study measured
survival, mortality, or viral load per animal, and a pathogen-challenge
contrast is not a tolerance phenotype. For how these feed candidate ranking, see
[implementation_status.md](implementation_status.md) and
[interpreting_candidate_scores.md](interpreting_candidate_scores.md).

## What this does not yet prove

The full-depth runs show the pipeline produces poolable, fully mapped
statistics from real reads, and that random-effects pooling works on real data.
They do not establish anything biological about these oysters, and they leave
plenty unverified:

* **Three comparisons are not run.** Of the three that are, Midori France has
  only 2 infected animals out of 5.
* **Two clean end-to-end runs.** The 2026-09-30 Midori Australia and Midori
  France runs are the first at full depth to finish `OK`, after `RENDER_REPORT`
  was moved to local scratch and the work directory off exFAT. Earlier runs
  all exited `ERR` at that step.
* **No tolerance phenotype.** Nothing here measures survival or mortality, so
  no evidence from this study can be `resilience_associated`, whatever the
  population contrast suggests. Viral load is now measured per animal, but
  from reads, as QC; it is not yet an analysed phenotype.
* **k = 2 is a thin pool.** Two studies give a very imprecise heterogeneity
  estimate (τ², I²). Candidates supported by this pool alone still need
  replication.
* **No like-for-like population contrast yet.** See
  [Why only one comparison pools](#why-only-one-comparison-pools).
* **No covariance-aware within-study model.** Using more than one comparison
  from this study in a single pool needs within-study covariance modelled;
  `meta_analysis_primary` avoids the question rather than answering it.
* **No containers.** All full-depth runs used native Homebrew tools. Every
  image tag in `containers/README.md` is still unverified, and reproducibility
  on another machine has not been demonstrated.
* **QC metrics are not in the manifests.** See [Full-depth runs](#full-depth-runs).
