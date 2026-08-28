# Candidate studies to curate

A screened shortlist of real, public *M. gigas* datasets to register after
`HESSER2024_VCOR`. Every accession here was verified against NCBI BioProject
and the ENA read-run API on **2026-08-28**; sample counts and group structures
come from the deposited run metadata, not from paper text.

This is a working backlog, not a commitment. Update it as studies are curated
or rejected.

## The selection criterion that actually matters

The obvious filter — "does the paper publish full statistics?" — is the wrong
one to lead with.

`HESSER2024_VCOR` harmonizes cleanly, resolves 87.2% of its identifiers, and
contributes **nothing** to any meta-analysis, because its supplementary table
reports only `log2FoldChange` and `padj`. No standard error, no unadjusted
p-value, so no inverse-variance weight. That is not an unlucky draw: a
significance-filtered table with adjusted p-values only is the *normal* format
for supplementary DE tables in this literature. Selecting studies by hoping
their supplementary files are richer will mostly reproduce the same dead end.

**The reliable route to poolable evidence is raw reads in SRA.** When AREE
reanalyzes raw data itself, its own DESeq2 run produces `lfcSE` and an
unadjusted p-value by construction — the study becomes poolable regardless of
what the authors chose to publish. It also finally exercises the
`raw_reanalysis` Nextflow path, which has never been run against real data
(see [implementation_status.md](implementation_status.md)).

So the screen applied here, in priority order:

1. **Raw reads deposited and public** (SRA/ENA), not just a supplementary table.
2. **Replicated design**, n ≥ 3 per group. See the rejected candidate below —
   this is not a formality.
3. **Resilience-relevant contrast.** Prefer a tolerance/resistance phenotype
   (resistant vs susceptible lineages) over pure exposure (treated vs control),
   per [resilience_vs_exposure.md](resilience_vs_exposure.md). Judge this on
   whether a phenotype was **measured**, not on how the title is worded — see
   the correction above.
4. **RNA-seq first**, because the intake converter currently handles only
   differential-expression tables and the RNA-seq workflow is the most complete.
5. **Stressor and phenotype spread**, so that pooled groups have more than one
   study in them. Two studies of the same stressor beat five of five different
   ones — a meta-analysis of k=1 is not a meta-analysis.

## Tier 1 — curate these first

> **#1 is registered** as `CALLA2026_OSHV` — design verified, sample sheet and
> FASTQ manifest generated, reanalysis not yet run. See
> [first_raw_reanalysis.md](first_raw_reanalysis.md).
>
> **Correction from the first draft of this page.** #1 was listed here as a
> study whose "contrast *is* a resilience phenotype, not an exposure", on the
> strength of its BioProject title (*"Evaluating Pacific oyster lineages for
> tolerance to Ostreid herpesvirus"*). Curating it showed that is not what the
> deposited data contains: there is no survival, mortality, or viral-load
> measurement for these animals, and the publication frames its results as
> groundwork for *future* tolerance breeding. It is registered as
> `disease_susceptibility`, not `disease_resistance`. A resilience-sounding
> project title is not a resilience phenotype — check for a measured outcome
> before promising one.

| # | BioProject | Runs | Design (from run metadata) | Resilience context | Why first |
|---|---|---|---|---|---|
| 1 ✅ | `PRJNA1329250` | 42 RNA-seq | 2 populations × 4 viral-strain levels (Control/Australia/France/USA), n=5 challenged / n=6 control | OsHV-1 challenge in two hatchery populations | Well replicated across all eight groups, and the same stressor class as the study already registered, so the two can eventually pool. Calla et al. 2026 ([10.1016/j.fsi.2026.111154](https://doi.org/10.1016/j.fsi.2026.111154)). |
| 2 | `PRJNA593309` | 43 RNA-seq | OsHV-1 × temperature (21/26/29 °C) × timepoint, n=3 | Disease resistance under thermal modulation | Multi-stressor, well replicated, and **published open access** — Delisle et al. 2020, *J Exp Biol* ([10.1242/jeb.226233](https://doi.org/10.1242/jeb.226233)). Pairs with #1 on pathogen challenge. |
| 3 | `PRJNA826964` | 18 RNA-seq | control vs OA × 3 timepoints (7/28/56 d), n=3 | Ocean acidification, energy metabolism | Clean 2×3 factorial, small enough to reanalyze quickly, and opens a second stressor class. |

Doing #1 and #2 together is the point: they give the pathogen-challenge group
**k ≥ 2 with real standard errors**, which is the first time random-effects
pooling would run on anything but simulated data.

## Tier 2 — good candidates, specific caveats

| # | BioProject | Runs | Design | Context | Caveat |
|---|---|---|---|---|---|
| 4 | `PRJNA678408` | 10 WGBS | 5 diploid + 5 triploid | Desiccation + acute heat | Would be AREE's **first real methylation study**. But the intake converter does not handle region tables yet, and coordinate-based evidence is exposed to the assembly change described in [handling_genome_versions.md](handling_genome_versions.md). |
| 5 | `PRJNA762441` | 18 RNA-seq | diploid/triploid × 3 timepoints, n=3 | Thermal stress, ploidy contrast | Clean design; ploidy is a useful resilience covariate. No linked publication confirmed. |
| 6 | `PRJNA913164` | 72 | diploid/triploid, marine heatwave | Thermal tolerance | **Tag-seq (3′ counts), not standard RNA-seq** — quantification differs from the salmon-based workflow. Largest n on the list. No publication found; treat as unpublished data. |
| 7 | `PRJNA735889` | 76 (RNA-seq + amplicon) | individually labelled, ~n=5/group | OA "tipping point" | Mixed assay project; the amplicon runs are a separate microbiome experiment and must be excluded at intake. |
| 8 | `PRJNA1196326` | 6 RNA-seq | 2 groups × n=3 | Transgenerational OA | Very small, but transgenerational designs are directly relevant to breeding and rare in the public record. |
| 9 | `PRJNA877226` | 6 RNA-seq | run labels blank in ENA | Vibrio × high temperature | Multi-stressor and cheap, but the deposited metadata does not describe groups — the design must be recovered from the paper before it is worth registering. |

## Adjacent, but a different evidence type

- `PRJNA1190893` — low-salinity adaptation, **102 WGS runs**. Population
  genomics, not differential expression. Registering it would mean a new
  harmonizer for variant/selection-scan evidence and a new `feature_type`.
  Worth doing eventually — salinity tolerance is a real breeding target and
  AREE currently has no real salinity evidence — but it is a feature, not a
  curation task.

- **Heat-resistant vs heat-susceptible families** (Baja California breeding
  program), Escobedo-Fregoso et al. 2023, *Comp Biochem Physiol D*
  ([10.1016/j.cbd.2023.101089](https://doi.org/10.1016/j.cbd.2023.101089)).
  The framing is exactly what AREE wants — RR vs SS phenotypes from a breeding
  program under oscillatory thermal challenge. I could not confirm a public
  raw-data accession for it; worth chasing, possibly by contacting the authors.

## Rejected, and why it matters

`PRJNA623063` — *"Transcriptome of the Pacific Oyster, Crassostrea gigas Larvae
after Vibrio alginolyticus Challenge"*, 12 RNA-seq runs.

On title alone this looked like the ideal second study: same life stage, same
stressor class, and nearly the same phenotype as `HESSER2024_VCOR`. The run
metadata says otherwise — all 12 runs share a single `sample_alias`
(`Oyster_M49`) across 12 timepoints (`T01`–`T12`). It is an **unreplicated time
course**, so no valid differential-expression contrast can be computed from it.

Recording this because the failure is invisible from the abstract, and because
a curator working from titles would have spent real effort before finding out.
Check `sample_alias` before downloading anything.

## Sourcing notes

- 93 candidate BioProjects were found for *C. gigas* / *M. gigas* across the
  stressor terms in AREE's ontology (thermal, acidification, pathogen,
  salinity, hypoxia, desiccation). The shortlist above is the subset passing
  the screen.
- GEO holds only 10 oyster expression series and is largely pre-2015; the
  useful corpus is in SRA/ENA. Do not treat a GEO search as sufficient.
- Hypoxia is badly underrepresented: exactly one BioProject matched, and it has
  no reads in ENA. AREE's phenotype ontology carries `hypoxia_tolerance`, but
  there may be no public transcriptomic dataset to populate it.

Bibliographic metadata for the two cited papers was retrieved from PubMed.
