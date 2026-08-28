# The first raw-reanalysis study

`CALLA2026_OSHV` (BioProject `PRJNA1329250`) is the second real study
registered in AREE and the first in `raw_reanalysis` mode. This page records
what was curated, the judgement calls made, and — plainly — what has not been
run.

**Status: registered, pipeline proven on a pilot, not harmonized.** The design,
sample sheet, FASTQ manifest, and reference context are complete and verified.
The RNA-seq workflow has now been executed end to end on a **subsampled** slice
of this study — see [What the pilot found](#what-the-pilot-found) — but no
full-depth run has been done, and **nothing from this study appears in the
evidence table or in any meta-analysis.**

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

# 3. Run one comparison. Keep -work-dir on APFS: Quarto's cleanup fails on exFAT.
nextflow run workflows/rnaseq -profile local \
  -c config/CALLA2026_OSHV.config \
  -work-dir /tmp/aree_nxf \
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

Then update `analysis_status` in the study YAML — `aree harmonize` warns if it
still reads `not_started`.

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

## What this does not yet prove

The pilot proves the pipeline runs and produces poolable statistics. It does
not prove anything biological about these oysters, and it leaves plenty
unverified:

* **No full-depth run.** One contrast at 1M read pairs per library is 2.7% of
  the data, for one of six comparisons.
* **No containers.** The run used native Homebrew tools. Every image tag in
  `containers/README.md` is still unverified, and reproducibility on another
  machine has not been demonstrated.
* **No evidence harmonized.** `CALLA2026_OSHV` still contributes zero records,
  and its `analysis_status` is still `not_started`.
* **Random-effects pooling on real data remains unexercised.** That needs a
  full-depth run here plus at least one more real study — `PRJNA593309` is the
  intended partner, see [candidate_studies.md](candidate_studies.md).
