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
2. **Replicated design**, n ≥ 3 per group. See the rejected candidates below —
   this is not a formality, and it is checked twice: distinct `sample_alias`
   values before download, then replicate dispersion after quantification
   (`scripts/check_replicate_dispersion.py`). Distinct sample names do not
   prove distinct animals.
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
> **#3 failed QC** as `IOCAS2022_OA_ENERGY` (registered 2026-09-30, reanalyzed
> 2026-10-01, QC decision 2026-10-07). Its replicates are technical, not
> biological — see *Rejected, and why it matters* below.
>
> **#7 is registered and harmonized** as `LUTIER2022_OA_TIPPING` (PR #11,
> merged 2026-10-09), AREE's first usable real ocean-acidification study.
> Two things on this page were wrong about it. First, the 31 "AMPLICON" runs
> are RNA-seq, mislabelled at deposit: all 76 runs make the paper's 15 tanks
> × 5 oysters. Second, the design is one tank per pH, so tank, not oyster, is
> the unit of replication. The prespecified contrast (pHT ≤ 6.7 vs ≥ 7.6,
> 4 vs 3 tanks) is fitted on pseudo-bulk tanks (`--replicate_unit tank`). It
> is classified `exposure_only`: the physiology was measured per tank, and no
> oyster was scored for tolerance.
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
| 3 ❌ | `PRJNA826964` | 18 RNA-seq | control vs OA × 3 timepoints (7/28/56 d), n=3 | Ocean acidification, energy metabolism | Clean 2×3 factorial, small enough to reanalyze quickly, and opens a second stressor class. |

Doing #1 and #2 together is the point: they give the pathogen-challenge group
**k ≥ 2 with real standard errors**, which is the first time random-effects
pooling would run on anything but simulated data.

## Tier 2 — good candidates, specific caveats

| # | BioProject | Runs | Design | Context | Caveat |
|---|---|---|---|---|---|
| 4 | `PRJNA678408` | 10 WGBS | 5 diploid + 5 triploid | Desiccation + acute heat | Would be AREE's **first real methylation study**. But the intake converter does not handle region tables yet, and coordinate-based evidence is exposed to the assembly change described in [handling_genome_versions.md](handling_genome_versions.md). |
| 5 | `PRJNA762441` | 18 RNA-seq | diploid/triploid × 3 timepoints, n=3 | Thermal stress, ploidy contrast | Clean design; ploidy is a useful resilience covariate. No linked publication confirmed. |
| 6 | `PRJNA913164` | 72 | diploid/triploid, marine heatwave | Thermal tolerance | **Tag-seq (3′ counts), not standard RNA-seq** — quantification differs from the salmon-based workflow. Largest n on the list. No publication found; treat as unpublished data. |
| 7 ✅ | `PRJNA735889` | 76 RNA-seq (31 mislabelled AMPLICON) | 15 tanks, one per pH, × 5 oysters | OA "tipping point" | Registered as `LUTIER2022_OA_TIPPING`. **Do not filter on ENA library_strategy**: the AMPLICON runs are RNA-seq. One tank per pH, so fit on tanks, not oysters. |
| 8 | `PRJNA1196326` | 6 RNA-seq | 2 groups × n=3 | Transgenerational OA | Very small, but transgenerational designs are directly relevant to breeding and rare in the public record. |
| 9 | `PRJNA877226` | 6 RNA-seq | run labels blank in ENA | Vibrio × high temperature | Multi-stressor and cheap, but the deposited metadata does not describe groups — the design must be recovered from the paper before it is worth registering. |

## Adjacent, but a different evidence type

- `PRJNA1190893` — low-salinity adaptation, **102 WGS runs**. Population
  genomics, not differential expression. Registering it would mean a new
  harmonizer for variant/selection-scan evidence and a new `feature_type`.
  Worth doing eventually — salinity tolerance is a real breeding target and
  AREE currently has no real salinity evidence — but it is a feature, not a
  curation task.

- **Heat-resistant vs heat-susceptible families** (CIBNOR breeding program,
  Baja California Sur), Arredondo-Espinoza et al. 2023, *Comp Biochem Physiol
  Part D Genomics Proteomics* 47:101089
  ([10.1016/j.cbd.2023.101089](https://doi.org/10.1016/j.cbd.2023.101089),
  PMID 37269757). Earlier drafts of this page cited it as "Escobedo-Fregoso et
  al." — she is the last and corresponding author; the first author is
  Arredondo-Espinoza. Author list and journal name verified against PubMed on
  2026-09-05.

  The framing is exactly what AREE wants — RR vs SS phenotypes from a breeding
  program under an oscillatory thermal challenge (26–34 °C, 30 days, sampled at
  day 0 and day 30). It would be AREE's **first real study with a measured
  resilience outcome**, and its first thermal stressor; every real candidate in
  the current pool is `disease_associated` for want of exactly this. Per-group
  replication is not stated in the abstract and must be confirmed before
  registration.

  No public raw-data accession could be located as of 2026-09-05: BioProject
  searches on organism plus thermal-challenge terms and on author/CIBNOR both
  returned nothing, and the Europe PMC record reports no associated data or
  supplement with no open-access full text to read the availability statement
  from. Note that **`sr320` is a co-author**, so lab storage and a direct ask to
  the corresponding author are better first moves than cold outreach. Tracked in
  [issue #6](https://github.com/sr320/AREE/issues/6).

  **Update 2026-09-30:** CIBNOR's WGBS project `PRJNA690951` covers the same
  RR/SS design (see M1 below), so the design has public raw data now. The
  RNA-seq reads for this paper are still not located.

## DNA methylation and proteomics candidates

Screened **2026-09-30**. Sources: ENA read-run search for every
`Bisulfite-Seq`/`MeDIP-Seq`/`MBD-Seq` run under taxon 29159 (20 BioProjects);
PRIDE Archive search for *C. gigas* / *M. gigas* (19 projects); PubMed title/abstract
search (37 methylation papers, 66 proteomics papers). Group sizes come from
run labels or BioSample attributes, not from paper text. A publication is
"confirmed" only where Europe PMC finds the accession in the paper's full
text, or PRIDE records it as the dataset reference. Otherwise the link is
inferred from the title and design, and must be checked at intake.

Neither assay can go straight into the processed-table converter. That
converter handles DE tables only (see
[implementation_status.md](implementation_status.md)). So methylation studies
need the raw Bismark/methylKit path, which has never been run on real data,
and proteomics studies need a per-study abundance-table script. Coordinate
evidence from methylation is also exposed to the assembly change in
[handling_genome_versions.md](handling_genome_versions.md).

### Methylation: Tier 1

> **M1 is registered** as `CIBNOR2021_HEAT_RRSS` (2026-09-30) — AREE's first
> real methylation study. All 12 libraries share one BioSample, so the design
> was parsed from `experiment_title` with the new `--label-field` /
> `--label-pattern` options of `aree fetch-samplesheet`. Reanalysis is **blocked**
> until the methylKit step models family: with 2 families per class, running it
> as-is would treat 6 vs 6 libraries as independent.

| # | Accession | Runs | Design (from deposited metadata) | Context | Why |
|---|---|---|---|---|---|
| M1 ✅ | `PRJNA690951` | 12 WGBS | 2 thermal-**resistant** (RR52, RR59) + 2 **susceptible** (SS05, SS35) families × n=3, gill, day 30 of 26–34 °C oscillation (CIBNOR) | Thermal tolerance, **measured phenotype** | Same breeding program and challenge as Arredondo-Espinoza et al. 2023, the RR/SS study in [issue #6](https://github.com/sr320/AREE/issues/6). This is the first public raw data from that design. No methylation paper cites it yet. It is the **methylation** arm: the RNA-seq reads are still not located. Effective replication is 2 vs 2 families, not 6 vs 6 oysters. |
| M2 | `PRJEB81880` | 40 EM-seq | 5 resistant vs 5 susceptible families, gill + mantle, before (T0) and after (T1) POMS, one oyster per family × tissue × time | Disease **resistance** (POMS) | A measured R/S phenotype with a before/after design, in AREE's best-populated stressor class (OsHV-1/POMS: `CALLA2026_OSHV`, `DELISLE2020_OSHV_TEMP`). Inferred to be Valdivieso et al. 2025, *Sci Total Environ* ([10.1016/j.scitotenv.2025.178385](https://doi.org/10.1016/j.scitotenv.2025.178385)). |
| M3 | `PRJNA562805` | 12 WGBS | intertidal vs subtidal origin × control/heat, n=3, gill | Thermal response by habitat origin | **Confirmed**: Wang et al. 2021, *Heredity* ([10.1038/s41437-020-0351-7](https://doi.org/10.1038/s41437-020-0351-7)). Pairs with M1 to give thermal methylation k=2. |
| M4 | `PRJNA682817` | 24 WGBS | diploid/triploid × pH 8.2/7.7, n=6 per cell, adult ctenidia (UW) | Ocean acidification (exposure) | Best-replicated OA methylation design found. Its natural RNA-seq partner, `IOCAS2022_OA_ENERGY`, failed QC; `LUTIER2022_OA_TIPPING` is OA but juvenile whole tissue at far lower pH, so not a like-for-like pairing. Lab-internal, so metadata recovery is easy. No publication located. Sister project of `PRJNA678408` (Tier 2 #4 above). |
| M5 | `PRJNA806944` | 8 WGBS | female gonad, ambient vs low pH | Ocean acidification (exposure) | **Confirmed**: Venkataraman et al. 2022, *BMC Genomics* ([10.1186/s12864-022-08781-5](https://doi.org/10.1186/s12864-022-08781-5)). Group labels are not in the run metadata, so recover the split from the paper. |

### Methylation: Tier 2

| Accession | Runs | Design | Caveat |
|---|---|---|---|
| `PRJEB105019` | 60 WGBS | 2 families (F14R, H2D) × 3 age cohorts at T0, plus post-challenge **R vs S** individuals, n=6 (Ifremer DECICOMP) | Best design on this list: individual-level resistance outcome, n=6. Released 2025-12 with no publication yet, so treat as unpublished and check the embargo/reuse terms. |
| `PRJNA609264` | 47 WGBS | early microbial exposure vs control, F1/F2, 2 families, n=3 at most timepoints | **Confirmed**: Fallet et al. 2022, *Microbiome*. Transgenerational disease protection is directly relevant to breeding. However, the design is a sprawling time course with pools and mixed assays, so intake needs a curated subset. |
| `PRJNA807732` | 24 WGBS | intertidal F0/F1/F2 and subtidal F0, control vs heat, n=3 | Probably Wang et al. 2023, *Sci Total Environ* (transgenerational intertidal). Not confirmed. |
| `PRJNA1113357` | 9 WGBS | *V. alginolyticus*, gill, 0/6/48 h, n=3 | Probably Li et al. 2024, *Fish Shellfish Immunol*. Small and exposure-only, but pairs with `HESSER2024_VCOR` on *Vibrio*. |
| `PRJEB58545` | 48 WGBS | pesticide mixture, F0 × F1 exposure (E/T), gastrula and metamorphosis, n=4 | Probably Sol Dourdin et al. 2024, *Environ Sci Technol*. Opens `pollutant_exposure`, but no resilience phenotype. |
| `PRJEB60400` | 246 WGBS | POMS-adapted vs naive populations | **Confirmed**: Gawra et al. 2023, *Sci Adv*. This is population-epigenetic differentiation, not a treatment contrast, so it needs the same new evidence type as the WGS salinity project above. |

Skipped: `PRJNA684407` and `PRJNA833956` (triploid infertility and growth,
with no stressor). Also skipped are `PRJNA324546` (developmental MeDIP) and the
2013–2016 reference methylomes (no contrast).

### Proteomics

The raw-data pool is thin. PRIDE holds 19 *C. gigas* projects, and most
Chinese-lab proteomics papers on the stressor list have no PRIDE deposit. They
may be in iProX, which was not searched. Harvesting processed tables from the
papers is the realistic route. The Roberts-lab projects make that easiest.

| # | Accession | Design | Context | Notes |
|---|---|---|---|---|
| P1 | `PXD002316` | 4 treatments (temperature × pH) × 3 tank replicates, larvae | OA × warming (multi-stressor) | Harney et al. 2016, *J Proteomics* (PRIDE reference). Complete submission. |
| P2 | `PXD015434` | early juveniles, temperature | Thermal | **Confirmed**: Crandall et al. 2022, *PeerJ*. |
| P3 | `PXD013262` | seed, time × temperature | Thermal | **Confirmed**: Wanamaker et al. 2020, *BMC Genomics*. Together, P2 and P3 give thermal proteomics k=2, and both pair with M1/M3 for multi-omics convergence. |
| P4 | `PXD000835` | ctenidia, ambient vs high pCO₂ | Ocean acidification | **Confirmed**: Timmins-Schiffman et al. 2014, *BMC Genomics*. An OA protein layer alongside `LUTIER2022_OA_TIPPING` (RNA-seq) and M4/M5. |
| P5 | `PXD011365` | gill mitochondria, hypoxia–reoxygenation | **Hypoxia tolerance** | Sokolov et al. 2019, *J Proteomics* (PRIDE reference). This is the only public omics dataset found for `hypoxia_tolerance`, which the RNA-seq screen flagged as empty. |
| P6 | `PXD065668` / `PXD065706` | gill proteome + phosphoproteome, ERK inhibition | Thermotolerance divergence | **Confirmed**: Wang et al. 2026, *Commun Biol*. This is a pharmacological perturbation, so it is mechanism support, not population evidence. |
| P7 | `PXD064651` | mucus, *Vibrio* infection, iTRAQ | Pathogen challenge | No publication. Partial submission. |
| P8 | `PXD000905` | low vs high OsHV-1 load, 2-DE | Pathogen | Corporeau et al. 2014, *J Proteomics* (PRIDE reference). This is 2-DE spot data, so feature mapping will be weak. |

No PRIDE deposit was found for these papers, so they are processed-table only.
Leprêtre et al. 2020, *Front Immunol* (OsHV-1 proteomics in two families with
**contrasting susceptibility**) is the strongest proteomic resilience phenotype
found. Kim et al. 2025, *CBP-D* covers hypoxia + heat in hemocytes, with
transcriptomics and proteomics together. Chen et al. 2023, *Fish Shellfish
Immunol* covers salinity tolerance with multi-omics.

**Rejected:** `PXD008057` (two populations × 0/6/24 h heat). Its sample
protocol pools 10 oysters per timepoint per population, which leaves six
pooled samples and no replicates. This is the proteomic version of
`PRJNA623063` below.

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

`PRJNA826964` — *"The compromised energy management of Pacific oysters
(Crassostrea gigas) under ocean acidification conditions"*, 18 RNA-seq runs.
Registered as `IOCAS2022_OA_ENERGY` and kept in the registry with
`qc_status: failed`, so the reason stays on record.

This one passed the `sample_alias` check: 18 distinct aliases, labelled
"biological replicate 1–3" in every group. It failed only after full-depth
reanalysis. Within each group, replicates vary about as little as read
sampling alone would make them — median dispersion 0.0002–0.0016, against
0.03–0.05 for `CALLA2026_OSHV`'s real biological replicates. Each group is
most likely one RNA sample sequenced three times. The symptoms, in the order
they appeared: 43–48% of genes "DE" at every timepoint; PC1 = 96.5% of
variance within a 3 vs 3 contrast; three OA groups that agree with the
controls as well as with each other; and fold changes that barely replicate
between timepoints.

The danger is not wasted compute. Harmonized, such a study carries an `lfcSE`
far too small, so inverse-variance pooling would let it outweigh every real
study it joins. Run `scripts/check_replicate_dispersion.py` on every
`raw_reanalysis` study before harmonizing it.

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

Bibliographic metadata for every cited paper was retrieved from PubMed — including the author order, which one of them had wrong until 2026-09-05.
