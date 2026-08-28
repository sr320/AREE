# The first raw-reanalysis study

`CALLA2026_OSHV` (BioProject `PRJNA1329250`) is the second real study
registered in AREE and the first in `raw_reanalysis` mode. This page records
what was curated, the judgement calls made, and — plainly — what has not been
run.

**Status: registered, not harmonized.** The design, sample sheet, FASTQ
manifest, and reference context are complete and verified. The pipeline run
that would produce result tables has not been executed. Nothing from this study
appears in the evidence table or in any meta-analysis.

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

Everything below is written and parameterized. None of it has been executed.

```bash
# 1. Stage the reads (226 GB), then verify against ENA's checksums.
awk -F'\t' 'NR>1 {print $4"  "$1"_"$2".fastq.gz"}' \
  data/studies/CALLA2026_OSHV/fastq_manifest.tsv > md5sums.txt
md5sum -c md5sums.txt

# 2. Fetch the reference transcriptome and annotation.
datasets download genome accession GCF_963853765.1 --include rna,gtf

# 3. Run one comparison.
nextflow run workflows/rnaseq \
  -config config/CALLA2026_OSHV.config \
  -profile docker \
  --comparison_id midori_oshv1_france_vs_control

# 4. Harmonize the standardized output it produces.
export AREE_CROSSWALK=data/reference/crosswalk/mgigas_gene_id_crosswalk.tsv
aree harmonize --study CALLA2026_OSHV \
  --input results/rnaseq/CALLA2026_OSHV/midori_oshv1_france_vs_control/standardized.tsv
```

Then update `analysis_status` in the study YAML — `aree harmonize` warns if it
still reads `not_started`.

## What this does not yet prove

The RNA-seq Nextflow workflow remains a **scaffold that has never been executed
against real FASTQ** ([implementation_status.md](implementation_status.md)).
Registering this study does not change that. The honest summary is that AREE
now has a real, verified, fully specified raw-reanalysis job queued, and the
first attempt to run it should be expected to surface bugs in the workflow —
that is the point of running it.

Until it runs, random-effects pooling on real data remains unexercised.
