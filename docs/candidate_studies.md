# Candidate studies to curate

A working backlog of real, public *M. gigas* datasets, initially screened
on **2026-08-28** and expanded in dated updates below. Sample counts and group
structures come from deposited metadata unless a paper is explicitly cited.
For current registration and harmonization status, use the [dashboard](index.qmd)
and study YAMLs; historical screening notes below are not a live status table.

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
what the authors chose to publish. The
`raw_reanalysis` Nextflow path has now run at full depth on real data
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

> **#1 is registered** as `CALLA2026_OSHV` — three of six comparisons have
> full-depth reanalysis and harmonized evidence. See
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
| 6 🟡 | `PRJNA913164` | 72 | diploid/triploid, marine heatwave | Thermal tolerance | **Tag-seq (3′ counts), not standard RNA-seq.** Run it with `--tag_seq` (no length correction in Salmon or DESeq2; see `workflows/rnaseq/README.md`). Largest n on the list. Registered as `GEORGE2023_HEATWAVE_PLOIDY` (#18): it **is** published, George et al. 2023, *Glob Chang Biol* ([10.1111/gcb.16880](https://doi.org/10.1111/gcb.16880)); the BioSamples carry no design, which comes from the authors' NCBI submission sheet. |
| 7 ✅ | `PRJNA735889` | 76 RNA-seq (31 mislabelled AMPLICON) | 15 tanks, one per pH, × 5 oysters | OA "tipping point" | Registered as `LUTIER2022_OA_TIPPING`. **Do not filter on ENA library_strategy**: the AMPLICON runs are RNA-seq. One tank per pH, so fit on tanks, not oysters. |
| 8 ❌ | `PRJNA1196326` | 6 RNA-seq | 2 groups × n=3 | Transgenerational OA | Registered as `IOCAS2024_OA_TRANSGEN` (#17) with `qc_status: failed`: each arm is one larval pool sequenced three times. See *Rejected* below. |
| 9 ❌ | `PRJNA877226` | 6 RNA-seq | infected only, 12 h vs 48 h, n=3 | Vibrio × high temperature | Rejected: no control arm. See *Rejected* below. |

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
need assay-specific raw processing; CIBNOR2021 has now completed Bismark plus
pooled-family limma analysis. Proteomics can use the quantified-matrix workflow,
as TIMMINS2014 does, with per-study intake for published tables. Coordinate
evidence from methylation is also exposed to the assembly change in
[handling_genome_versions.md](handling_genome_versions.md).

### Methylation: Tier 1

> **M1 is registered** as `CIBNOR2021_HEAT_RRSS` (2026-09-30) — AREE's first
> real methylation study. All 12 libraries share one BioSample, so the design
> was parsed from `experiment_title` with the new `--label-field` /
> `--label-pattern` options of `aree fetch-samplesheet`. Reanalysis completed
> on 2026-10-16 using family-level pooling and limma: 2 vs 2 families, not
> 6 vs 6 independent oysters. The 132,587 harmonized regions yielded no
> BH-significant differences; this is a low-power null result.

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
| `PRJEB105019` | 60 WGBS | 2 families (F14R, H2D) × 3 age cohorts at T0, plus post-challenge **R vs S** individuals, n=6 (Ifremer DECICOMP) | Best design on this list: individual-level resistance outcome, n=6. Released 2025-12 with no publication yet, so treat as unpublished and check the embargo/reuse terms. The EM-seq arm of DECICOMP: its metabolomics arm is `MTBLS13463`, registered as `DECICOMP2025_POMS_AGE` (#21). Run titles give the cohort codes: C3 = 4, C7 = 16, C8 = 28 months. Its 24 resistant/susceptible libraries are all from the 28-month cohort, at 500 hpc. |
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
| P1 | `PXD002316` | 4 treatments (temperature × pH) × 3 tank replicates, larvae | OA × warming (multi-stressor) | **Registered, not harmonized** as `HARNEY2016_OA_WARM_LARVAE` (2026-10-17). The design is sound, but this is a 2-DE study with no deposited spot volumes or gel images; PRIDE holds only spot identifications. The paper reports 12 of ~650 spots, selected at P < 0.02, with no SE, and 9 of 11 accessions are now deleted from UniProtKB. Worth reanalyzing if the SameSpots data are ever released. |
| P2 | `PXD015434` | early juveniles, temperature | Thermal | **Registered, QC failed** as `CRANDALL2022_HEAT_SEED` (2026-10-17): one silo per temperature, and the two silos differed beyond temperature (see Rejected below). |
| P3 | `PXD013262` | seed, time × temperature | Thermal | **Registered, QC failed** as `WANAMAKER2020_HEAT_SEED` (2026-10-17): one silo per temperature, technical replicates only (see Rejected below). With P2 also failed, no public thermal proteomics study in this list is analyzable. |
| P4 | `PXD000835` | ctenidia, ambient vs high pCO₂ | Ocean acidification | **Reanalyzed** as `TIMMINS2014_OA_GILL` (2026-10-17), the first real proteomics study in AREE. 4 individual oysters per pCO₂ × mechanical-stress group, from the paper's per-oyster NSAF table. Low-power null: 639 proteins testable, none at padj < 0.05. It gives an OA protein layer with SEs alongside `LUTIER2022_OA_TIPPING` (RNA-seq) and M4/M5. QC passed with warnings. |
| P5 | `PXD011365` | gill mitochondria, hypoxia–reoxygenation | **Hypoxia tolerance** | **Registered, not reanalyzed** as `SOKOLOV2019_HYPOX_MITO` (2026-10-17). It is the only public omics dataset found for `hypoxia_tolerance`, but not usable as deposited: the phosphoproteome has no run-to-condition key; the total proteome is 3 oysters per condition, acquired in condition blocks; and the only quantitative files are Scaffold mzIdentML exports over oyster_v9 RefSeq IDs. Reopen if the authors supply the phospho key and per-oyster NSAF tables. |
| P6 | `PXD065668` / `PXD065706` | gill proteome + phosphoproteome, ERK inhibition | Thermotolerance divergence | **Registered, not reanalyzed** as `WANG2026_ERK_HEAT` (2026-10-17). The quantification is only in proprietary Spectronaut `.sne` archives, and the paper's per-pool tables list only significant features. Proteins are the authors' own gene models, and the runs were acquired in condition blocks. Reopening needs a DIA re-search of the `.raw` files. It is a pharmacological perturbation, so mechanism support, not population evidence. |
| P7 | `PXD064651` | mucus, *Vibrio* infection, iTRAQ | Pathogen challenge | **Registered, cannot be reanalyzed** as `JIN2025_VIBRIO_MUCUS` (2026-10-17). IBT 16-plex with no channel-to-sample key, no quantification table, and a Mascot search run without quantitation. The infected and control channels cannot be told apart. Unpublished. |
| P8 | `PXD000905` | low vs high OsHV-1 load, 2-DE | Pathogen | **Registered, not harmonized** as `CORPOREAU2014_OSHV_LOAD` (2026-10-17). It is in the same situation as P1. PRIDE holds only identifications of the 25 excised spots; no spot volumes or gel images were deposited. The paper reports only average fold changes for those significant spots. |

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

**Registered, QC failed:** `PXD013262` (P3, `WANAMAKER2020_HEAT_SEED`).
The 26 runs are 13 pools in technical duplicate, and the authors' sample
metadata on figshare puts every 23 °C pool in one silo and every 29 °C pool in
another. Temperature is confounded with silo, and no replicate is biological.
The paper's ASCA found no temperature effect (P = 0.19). The authors' code
also mentions a second 23 °C silo that was never deposited. If those samples
exist, re-screen.

**Registered, QC failed:** `PXD015434` (P2, `CRANDALL2022_HEAT_SEED`).
The 16 runs are 4 pools × 4 DIA m/z-range injections. The pools come from
two silos, the only 23 °C control and the 29 °C control, each sampled on two
days. The hatchery notebook in the authors' Zenodo release shows the 23 °C silo
started at 16–18 °C, set poorly and died out after a ciliate outbreak. The
paper's survival contrast (0 % vs 86 %) is therefore confounded with that
history, and no protein reaches adjusted P < 0.05. Both P2 and P3 come from
the same Taylor Shellfish 23-vs-29 °C silo trials. Any further proteomics from
that program is worth checking for silo replication before claiming it.

2-DE studies (P1, P8) are a separate problem. PRIDE deposits for gel studies
usually hold only identifications of excised spots, not spot volumes, so they
cannot be requantified. Check for deposited gel images or a spot-volume table
before claiming one.

## Non-coding RNA, chromatin, and metabolomics candidates

Screened **2026-10-11**. Sources:
- every ENA read run under taxon 29159 (12,994 runs), bucketed by `library_strategy`;
- every `RNA-Seq` run with `library_selection = size fractionation`, or with miRNA/lncRNA/circRNA in its study or experiment title;
- MetaboLights (EBI search for "gigas");
- Metabolomics Workbench (species and title search);
- a PRIDE refresh;
- PubMed searches for *C. gigas* lncRNA/circRNA and miRNA under stress (18 papers).

Group sizes come from `sample_alias` and `experiment_title`.

Non-RNA-seq run counts for *M. gigas* in ENA on this date:

| Strategy | Runs |
|---|---|
| Bisulfite-Seq | 607 |
| miRNA-Seq | 117 |
| RIP-Seq | 44 |
| ChIP-Seq | 28 |
| ATAC-seq | 24 |
| MeDIP-Seq | 21 |
| ncRNA-Seq | 8 |

`AMPLICON` (3,329) and `OTHER` (2,069) runs are mostly microbiome and
genotyping runs, but `LUTIER2022_OA_TIPPING` shows that mislabelled RNA-seq can
hide there.

The methylation and proteomics findings above still hold: no new
stress-relevant WGBS/EM-seq project has appeared since 2026-09-30, and the one
new *C. gigas* PRIDE project (`PXD082789`, simulated digestion) is food science.

### Metabolomics: the real find of this screen

AREE's only metabolomics study is the simulated `GIGAS_GROW06`. Both of these
are POMS studies, AREE's best-populated stressor class, and both have a
**measured disease outcome**:

> **Status 2026-10-11.** X1 and X2 are registered as `DECICOMP2025_POMS_AGE`
> (#21) and `DUPERRET2025_POMS_STARVE` (#22). Both are draft PRs, because
> **neither can be analyzed**: every Metabolite Assignment File (MAF) in both
> deposits lists metabolite names but holds **no abundance values**. The only
> data are Agilent targeted-dMRM vendor files (`.d`), which no AREE workflow
> reads. Ask the submitters (J. Vidal-Dupiol, IHPE) to complete the MAFs, or
> add an MRM integration step. Curation also found, in the order it matters:
>
> - **X2: diet is fully confounded with MS batch** on all four platforms.
>   Fed samples were acquired Dec 2023–Feb 2024 and starved samples May–Jun
>   2024, with no overlap. Even with abundances, a diet effect here cannot be
>   separated from batch.
> - **X2 samples are pools of 10 oysters per tank × family**, with 3 tanks per
>   diet, so replication is 3 tanks.
> - **X1 age is confounded with challenge year**: each age cohort was
>   challenged in a different year. MS batch is not confounded with age.
> - **X1 sampling time (0 hpc) is inferred**, from the EM-seq arm's identical
>   2 × 3 × 6 structure at 0 hpc. Whether the metabolomics and EM-seq samples
>   are the **same animals is not stated**. The first draft of this table said
>   they were; that was an assumption.
> - Neither deposit has database identifiers or MSI levels, so metabolite
>   mapping confidence starts at `unresolved`.
> - The studies may not be independent: both use family F14R from the same
>   2020 Ifremer production, and X1's MAFs carry X2's sample columns.

| # | Accession | Design | Context | Notes |
|---|---|---|---|---|
| X1 🟡 | `MTBLS13463` (+ `PRJEB105019`) | DECICOMP (Ifremer): 2 families × 3 ages (4/16/28 mo) × 6 whole animals, 0 hpc before OsHV-1 cohabitation; older oysters are **more tolerant** of POMS | Disease **resistance**, measured per cohort | Registered as `DECICOMP2025_POMS_AGE` (#21): 28 vs 4 mo (primary) and 16 vs 4 mo. Same project as `PRJEB105019` (methylation Tier 2), which is cited but kept out of `accessions` so it can be registered on its own. Manuscript in preparation; no DOI. |
| X2 🟡 | `MTBLS13313` | 2 families × fed/starved (3 tanks each) × 0–12 hpc, pools of 10; starvation **reduced POMS mortality** (F14R 34.8% vs 75.5%; F11N 93.3% vs 100%) | Disease resistance × `nutritional_limitation` | Registered as `DUPERRET2025_POMS_STARVE` (#22): starved vs fed at 0 hpc. Published: Duperret et al. 2025, *iScience* 29:114333 ([10.1016/j.isci.2025.114333](https://doi.org/10.1016/j.isci.2025.114333)). Its RNA-seq arm is public and is the better study; see below. |

Beyond the empty MAFs, there is no converter for real MTBLS `m_*.tsv`
tables: metabolomics harmonization exists only for the demo format.

### The RNA-seq arm of X2

The "integrative omics" in X2's paper includes RNA-seq, deposited as two ENA
projects split by diet: `PRJEB86618` (starved) and `PRJEB86525` (fed). The
design is 2 families × 2 diets × 7 times (0–72 hpc) × 3 tank pools, 83 of 84
libraries deposited, single-end 75 bp. Unlike the metabolomics, library
numbering interleaves the two diets.

> **Registered** as `DUPERRET2025_STARVE_RNASEQ` (#23), one study across both
> projects. All three comparisons now have harmonized evidence. Comparisons:
> starved vs fed at 0 hpc (primary, `resilience`), matching X2's metabolomic
> contrast, and 24 vs 0 hpc within each diet (`disease`). It needed two
> workflow changes, now merged: single-end input (#25) and, as a general
> feature, `--covariates` (#27).
>
> **It is fitted on tanks, not with a family covariate.** Both families share
> each tank, so `~ family + condition` would count 6 libraries per diet when 3
> tanks were treated. Each tank's two family libraries are summed
> (`--replicate_unit tank`, or `tank_hpc` for the time contrasts), giving
> 3 vs 3. This is the same lesson as `LUTIER2022_OA_TIPPING`: find the unit
> the treatment was applied to before counting replicates.
>
> It shares tanks with X2, so the two are convergent evidence across omics
> layers, not independent replication.

Metabolomics Workbench has no *C. gigas* study. Its oyster entries are
*C. virginica* predator/shell studies.

### Chromatin accessibility (ATAC-seq)

| # | Accession | Runs | Design | Context | Notes |
|---|---|---|---|---|---|
| C1 | `PRJNA970954` | 6 ATAC + 6 RNA-seq | adult gill, heat (HC) vs no heat (NC), n=3, sample labels shared across both assays (HC.gigas1–3 / NC.gigas1–3) | Thermal (exposure) | Two-species study (*C. gigas* vs *C. angulata*). Paired ATAC + RNA-seq gives a clean, linkable multi-omics pair (peak → nearest gene → DE). Needs a new `feature_type` for accessible regions and a peak-calling workflow. |
| C2 | `PRJNA1212587` | 6 ATAC | intertidal vs subtidal, n=3 | Thermal / habitat origin | Probably the same group as M3 (`PRJNA562805`, Wang et al. 2021, intertidal vs subtidal WGBS). Would give methylation + accessibility on one contrast. No publication confirmed. |
| — | `PRJNA1113803` | 4 ATAC | fast vs slow growth strains, **n=2** | Growth | Skip: n=2, no stressor. |
| — | `PRJNA417263` | ATAC + 1 RNA | hepatopancreas, 2 technical reps | — | Skip: reference chromatin map, no contrast. |

### ChIP-seq / CUT&Tag: annotation, not evidence

`PRJNA577897` (HSF1 ChIP, gill, heat), `PRJNA1292992` (BTF3L4 CUT&Tag, HSP
regulation, 2026), and `PRJNA1193684` (CLOCK DAP-seq, thermal adaptation) are
all **n=1 vs input**. They cannot produce effect sizes. They can, however,
annotate candidates: "HSF1 binds this promoter" is useful mechanism context on
an evidence card. Treat them as a reference layer (`data/reference/`), not as
studies. `PRJNA1292166` is a 2026 ChIP method-optimization set (H3K4me3, gill
and mantle) and has the same use.

### Non-coding RNA

The deposited small-RNA record is weak for AREE's purposes. Most miRNA-Seq
projects are atlases, developmental series, or reproduction/triploidy studies.
The stress-relevant ones are mostly unreplicated:

| Accession | Runs | Design | Verdict |
|---|---|---|---|
| `PRJNA836789` | 6 miRNA + 6 ssRNA | selectively bred strains HH vs CC, n=3, soft tissue | **Out of scope (checked 2026-10-15).** The BioProject describes it as a study of **inbreeding depression** (one inbred strain with growth depression), not a stressor or resilience contrast. |
| `PRJNA750516` | 2 miRNA + 6 lncRNA-labelled RNA | two groups (RZ_GH vs KTD_GL), **one sample alias each** | Skip: replicates are re-sequenced libraries of one sample. |
| `PRJNA694837` | 6 miRNA + 6 RNA | BM vs OM mantle, n=3 | Probably shell colour/biomineralization. No stressor, so low priority. |
| `PRJNA341808` ❌ | 3 "miRNA" + 9 RNA | larvae: control / bacteria / heat / bacteria+heat, Ion Torrent | **Rejected 2026-10-14** (see *Rejected* below): the replicates are technical, and the "miRNA-Seq" control runs are mislabelled RNA-seq. |
| `PRJNA484109` | 1 miRNA + 9 RNA | OsHV-1 infection (Rosani et al. 2020, *BMC Genomics*, [10.1186/s12864-020-07026-7](https://doi.org/10.1186/s12864-020-07026-7)) | Only one small-RNA run is in this project. The paper's paired miRNA/mRNA libraries must be located before this is usable. Checked 2026-10-15: the six gill RNA-seq runs (`S1`–`S6`) carry no group labels in the deposit, so the design is not recoverable from it. |
| `PRJNA266361`, `PRJNA259995`, `PRJNA232734` | 1–3 | miRNA stress/atlas, single pooled sample | Skip: unreplicated. |

Several published stress miRNA studies have **no accession in their
abstracts**. Check their data-availability statements before treating them
as raw-data candidates:
- exosomal miRNAs after heat or *Vibrio* (2024, *Dev Comp Immunol*, [10.1016/j.dci.2024.105174](https://doi.org/10.1016/j.dci.2024.105174));
- miRNAs in haemocyte immune priming (2026, *Front Immunol*, [10.3389/fimmu.2026.1753252](https://doi.org/10.3389/fimmu.2026.1753252));
- haemocyte miRNAs after bacteria and heat (2014, *PLoS One*, [10.1371/journal.pone.0088397](https://doi.org/10.1371/journal.pone.0088397)).

**lncRNA and circRNA almost never have their own deposits.** The *C. gigas*
lncRNA stress papers found re-quantify ordinary RNA-seq:
- heat, built from 27 public RNA-seq sets ([10.1007/s10126-022-10140-0](https://doi.org/10.1007/s10126-022-10140-0));
- OsHV-1 antiviral response ([10.1016/j.fsi.2018.07.032](https://doi.org/10.1016/j.fsi.2018.07.032));
- growth at low temperature ([10.1186/s12864-023-09543-7](https://doi.org/10.1186/s12864-023-09543-7)).

The realistic route to lncRNA evidence is therefore a **workflow feature, not
a curation task**: add lncRNA loci to the quantification index and run them
through the existing RNA-seq path for studies already registered
(`CALLA2026_OSHV`, `DELISLE2020_OSHV_TEMP`, `LUTIER2022_OA_TIPPING`).
PolyA-selected libraries will under-sample non-polyadenylated lncRNAs and miss
circRNAs entirely, so record library selection as a caveat.

### Other non-RNA-seq projects noted, not shortlisted

- `PRJNA894789` — 5-Aza-dC demethylation × heat, gill, n=3 (0 h control /
  6 h H₂O / 6 h 5-Aza). The runs are RNA-seq despite the `size fractionation`
  label. It is a pharmacological perturbation of methylation, so it is
  mechanism support for M1/M3 rather than population evidence, like P6.
- `PRJNA291364` — diuron-exposed parents → spat, WGBS + RNA-seq, n=2. This is
  transgenerational `pollutant_exposure`, and too small to pool.
- `PRJNA689936` — "wild oysters" WGBS, 6 runs from **one** sample alias. Skip.
- `PRJNA748038` — m6A RIP-seq (MeRIP) across development. It has no stressor,
  and RNA modification is a new evidence type.
- `PRJEB77044` — 96 runs, BC field stations × mantle/hemolymph. Labelled
  `RNA-Seq` with `size fractionation`; whether it is small RNA or standard
  RNA-seq was not checked. If standard, it belongs on the RNA-seq list as a
  field/environmental gradient design.

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

`PRJNA1196326` — *"Transcriptional acclimation following transgenerational
exposure of oysters to ocean acidification"*, 6 RNA-seq runs (IOCAS).
Registered as `IOCAS2024_OA_TRANSGEN` (#17) with `qc_status: failed`. It is
the same failure as `PRJNA826964`, from the same institute: dispersion
0.0008–0.0010 after a full-depth run, 42.6% of genes "DE", and an
expressed-SNP screen showing each arm is **one larval pool** sequenced three
times (within-arm allele-frequency difference 0.061 against 0.051 for one
RNA). The two arms also descend from different parents, so parental pH is
confounded with parental genotype.

`PRJNA341808`: *C. gigas* D-hinged larvae, control / bacteria / heat /
bacteria + heat, 3 Ion Torrent runs per group, 26.2 GB. **Rejected before
download, 2026-10-14**, by `scripts/prescreen_replicates.py`. It is the
first study turned away by the pre-screen, not after a full reanalysis.

Each group's three runs share one BioSample (`SAMN05728707`–`10`) and one
`sample_alias`, which was already a warning. The pre-screen, run on the first
1 M reads per run from NCBI, gave these median dispersions:

| Group | Median dispersion |
|---|---|
| heat | −0.0010 |
| bacteria | 0.0039 |
| control | 0.0050 |
| bacteria + heat | 0.0058 |

Two groups fail and two sit at the 0.005 threshold. The pre-screen reads high
(LI2022_HEAT_PLOIDY's real replicates scored 0.019–0.022 on it), so all four
groups are in the technical-replicate range. Each group is most likely one
larval RNA sample sequenced three times.

Two things the screen also settled:
- The three control runs deposited as `miRNA-Seq` are ordinary RNA-seq. They
  map to the transcriptome as well as the stress runs (55% vs 53–58%).
- A genotype check could not have decided this case. Larval pools from one
  cross share parental genotypes, whether or not the pools are independent.

An expressed-SNP screen run independently on 2026-10-13 (first ~2 M reads
of the control and heat runs, 167 sites) agrees but cannot decide alone:
runs within a group differ in allele frequency by 0.062 against 0.057
expected from one RNA sample, and groups differ from each other by 0.146.
That is what one RNA per group looks like, but large pools from one cross
would look much the same, as noted above.

`PRJNA694496` — Tan et al. 2021, *Front Physiol*
([10.3389/fphys.2021.663023](https://doi.org/10.3389/fphys.2021.663023)),
larvae of a thermotolerance-selected line vs a control line, 6 runs.
Rejected because the design **cannot be recovered**. All six runs share one
BioSample, whose tissue ("gill") also contradicts the experiment title
("larva"). The run aliases are vendor IDs (`W056-01-T01`…`T06`), and the
paper says only "three replicates for each group" without mapping runs to
groups. Assigning groups would be a guess. The authors also note that their
survival difference "may not fully reflect the thermotolerance" of the lines.

`PRJNA877226` — *"RNA-seq of Pacific oyster infected with Vibrio pathogen
under high temperature"*, 6 runs (Ocean University of China). The groups are
recoverable after all: the BioSample titles give *V. alginolyticus*-infected
gill at high temperature, **12 h vs 48 h**, n=3. But the deposit has **no
uninfected control and no normal-temperature arm**. The run names (`C33-n`)
match the 23 °C arm (`C23-n`) of the same lab's `PRJNA756710`
(`LI2022_LOWSAL_VIBRIO`, #26). Using that study's 30 ppt libraries as a
control would cross BioProjects and probably sequencing batches, so it is not
done. A within-infection time contrast is not stress-vs-control evidence.

`PRJNA298285`: larval OA × warming, University of Brest (2015). It has 24 RNA-seq
runs: trochophore / early / late veliger × pH 7.9 / 8.1 × 20 / 22 °C, two
runs per group. Each group's two runs share **one BioSample**. **Rejected
before download, 2026-10-15**, by `scripts/prescreen_replicates.py` (first
1 M reads per run from NCBI, the four 24 h groups). Every group's median
dispersion is **negative** (−0.0047 to −0.0065): the two runs vary *less*
than read sampling alone predicts, which is the signature of one library
sequenced on two lanes. Effective n is 1 per group, for all 12 groups.

**Two lessons from these rejections.** One BioSample per group, or a whole
study under one BioSample, is the cheapest early warning of
pseudo-replication. And replication can be checked before the full download:
`scripts/prescreen_replicates.py` does it from a subsample. An expressed-SNP
genotype check (`check_replicate_genotypes.py`, currently in sr320/AREE-co)
complements it when libraries are single animals or small pools of unrelated
animals. It cannot decide for large pools of siblings, which converge on
their parents' allele frequencies whether or not they are independent.

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
