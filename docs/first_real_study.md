# Worked example: registering and harmonizing a real study

Every other page in this documentation describes the design. This one records what
actually happened the first time AREE was pointed at a real published study —
including the things that broke. It is the most useful page to read before
curating your own.

## The study

**Hesser J, Mueller RS, Langdon C, Schubiger CB (2024).** *Immunomodulatory
effects of a probiotic combination treatment to improve the survival of Pacific
oyster (Crassostrea gigas) larvae against infection by Vibrio coralliilyticus.*
Frontiers in Immunology 15:1380089.
[doi:10.3389/fimmu.2024.1380089](https://doi.org/10.3389/fimmu.2024.1380089)
(open access, CC BY 4.0).

Registered as [`HESSER2024_VCOR`](../registry/studies/HESSER2024_VCOR.yaml).

It was chosen because it is a genuine resilience study — larval survival under a
pathogen challenge, with and without a protective intervention — and because its
data are actually obtainable.

## Finding a usable study is the hard part

A survey of every *C. gigas* series in GEO (31 series at time of writing) found
**not one** that ships a differential-expression table. GEO holds raw reads,
microarray `.gpr` files, count matrices, and bedgraphs. NCBI-generated count
matrices do not cover this organism. The count matrices that do exist mostly use
legacy or assembly-specific identifiers (`Cg19979`, `XLOC_036894`, `CGI_10001851`)
rather than NCBI GeneIDs.

The practical consequence: for *C. gigas*, **processed-results harmonization
usually means a supplementary table from an open-access paper**, not a
repository download. Plan curation around that.

This study deposited no raw reads at all — its only public data are the
supplementary DEG tables — which is exactly the case
`processed_results_harmonization` exists for.

## Getting from the publisher's file to an AREE result file

The published artifact is an Excel workbook with three sheets, one per
comparison. `src/intake/supplementary_table.py` performs the minimum mechanical
reshaping and records what it did in
`data/studies/HESSER2024_VCOR/intake_provenance.json`, with checksums of both the
source file and every generated result file. It deliberately does not alter
identifiers, impute missing statistics, or filter beyond what the authors did.

Two things worth knowing:

* **The published sheet contains a repeated header row part-way down** (source
  row 71 of `18hrPB+Vc v LOnly`). Real supplementary tables are not clean tables.
  The converter coerces the numeric columns and drops non-numeric rows, counting
  them in the provenance record rather than letting the string `"log2FoldChange"`
  reach the harmonizer as an effect size.
* **One of the three published comparisons was deliberately not registered.**
  The "18 hr PB Only" arm (probiotic, no pathogen) involves no stressor, and the
  stressor ontology has no term for a beneficial-microbe exposure. Forcing it into
  `pathogen_challenge` would have been false. This is recorded in the study's
  `limitations` rather than quietly dropped — a gap in the vocabulary is a finding,
  not something to paper over.

## What the crosswalk actually did

351 evidence records, resolved against
`data/reference/crosswalk/mgigas_gene_id_crosswalk.tsv`:

| mapping_confidence | records | what these are |
|---|---:|---|
| `exact` | 274 | direct hits on the current NCBI annotation |
| `inferred` | 32 | GeneIDs NCBI has retired, remapped to their current replacement |
| `unresolved` | 45 | genes NCBI discontinued **with no replacement** |

**87.2% resolved.** Every single unresolved identifier was checked against NCBI
`gene_history`: all 35 distinct unresolved LOC identifiers are genes NCBI itself
withdrew without a successor. Nothing resolvable was left unresolved.

The 32 `inferred` records exist only because of the retired-GeneID sidecar. The
authors aligned to assembly `GCA_902806645.1` (cgigas_uk_roslin_v1) and its
2020-era annotation, while the crosswalk is built on the 2024 annotation; without
retired-ID remapping those 32 records would have been silently lost. This is the
annotation-version problem in [handling_genome_versions.md](handling_genome_versions.md)
showing up on the first real study rather than in theory.

The one remaining unresolved identifier is not a LOC ID at all: the published
table mixes identifier systems and lists a bare gene symbol, `Myd88`. It is left
`unresolved` on purpose. *M. gigas* has at least five MyD88-family genes
(`Myd88`, `Myd88-2`, and several unnamed paralogs), so automatically assigning a
bare symbol to one of them would manufacture a paralog-level claim the source
does not support. Cases like this belong in
`data/mappings/ambiguous_symbol_map.yaml` as a curated, documented decision.

## Three bugs the synthetic data had hidden

Real data broke things that 37 passing tests had not:

1. **90% of the study failed to resolve.** The RNA-seq harmonizer labels a
   `gene_id` column as `ncbi_gene_id`, and the crosswalk's `ncbi_gene_id` column
   holds *numeric* GeneIDs — but published tables put the `LOC` form there. Only
   retired IDs resolved, because that lookup happened to strip the prefix.
   `ncbi_gene_id` and `locus_id` are now searched interchangeably, which is
   correct by NCBI's own convention.
2. **Decorated identifiers.** The source reports
   `gene-LOC105320749|LOC105320749` — a GFF feature prefix plus a pipe-joined
   symbol. `identifier_candidates()` now strips this deterministic decoration for
   *lookup only*; `feature_id_original` keeps the published string verbatim.
3. **The manifest lied about provenance.** It recorded a hardcoded path to the
   demo crosswalk regardless of which crosswalk was used. Manifests now record
   the crosswalk actually in force, with its SHA-256 and that of the retired-ID
   table.

## Keeping real and simulated evidence apart

Once a real study exists, the demo studies are no longer harmless. Two changes
enforce the separation:

* **`simulated` is now a column in the evidence schema and part of the
  meta-analysis grouping key**, so a pooled estimate can never combine fabricated
  demo evidence with evidence from a real study. The Streamlit interface has a
  matching "Data origin" filter.
* **Each study chooses its own crosswalk from its `simulated` flag**
  (`crosswalk_for_study()`). A simulated study always resolves against the demo
  crosswalk; a real study refuses to run at all unless a real crosswalk is
  selected. This matters because the demo LOC numbers collide with real GeneIDs
  denoting different genes — harmonizing the demo studies against the real
  crosswalk attaches real NCBI GeneIDs to fabricated effect sizes.

## The study contributes no pooled estimate — by design

`aree meta-analyze` returns 34 pooled rows, none of them from this study.

The source publishes a log2 fold change and a Benjamini-Hochberg adjusted
p-value, and **no standard error and no unadjusted p-value**. There is no
defensible way to recover a standard error from an adjusted p-value, so AREE
declines to pool rather than imputing one. The evidence records are still in the
table, still searchable, still exportable, and still available to any future
analysis that does not need inverse-variance weighting.

This is the correct outcome, and it is worth stating plainly: **a real study can
be fully harmonized and still contribute nothing to a meta-analysis.** The
limitation was predicted in the study's `limitations` field before the pipeline
was run.

## Reproducing it

```bash
aree build-crosswalk                       # once; ~380 MB streamed from NCBI
export AREE_CROSSWALK=data/reference/crosswalk/mgigas_gene_id_crosswalk.tsv

# The published spreadsheet -> the two result tables AREE harmonizes.
aree intake-supplementary data/studies/HESSER2024_VCOR/intake.yaml --check

aree validate-study registry/studies/HESSER2024_VCOR.yaml
aree register-study  registry/studies/HESSER2024_VCOR.yaml --update
aree harmonize --study HESSER2024_VCOR
```

The result files are committed, so the harmonize step runs without re-downloading
anything from the publisher. The published source spreadsheet is committed too,
which is what lets the intake step be verified rather than trusted: `--check`
regenerates the tables into a temporary directory and compares checksums against
both the committed files and `intake_provenance.json`. Drop `--check` to rewrite
them; a re-run producing byte-identical output preserves the original
`date_generated`, so verifying reproducibility never shows up as a diff.

Reading `.xls`/`.xlsx` sources needs the optional spreadsheet engines:
`pip install -e ".[intake]"`. CSV and TSV sources need nothing extra.

### What the intake step deliberately will not do

The conversion is mechanical: select and rename columns, drop rows with no
identifier or no effect size, drop non-numeric rows (the published sheet
`18hrPB+Vc v LOnly` contains a repeated header part-way down), and record all
of it. It never imputes a missing statistic. Because this source reports only
`log2FoldChange` and `padj`, the `lfcSE` and `pvalue` columns come out **empty**,
and that emptiness is what later stops the study from being pooled. Turning
that gap into a plausible-looking number in a spreadsheet would have hidden the
single most important fact about this study.

Two guards keep that honest: the source `sha256` in the intake config is
verified before every conversion, and `scripts/check_real_study_evidence.py`
(run by CI) fails if `standard_error` or `p_value` ever becomes populated for
this study.
