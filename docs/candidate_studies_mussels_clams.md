# Candidate mussel and clam studies to curate

A screened shortlist of real, public mussel and clam RNA-seq datasets, the
backlog for WP5 of [mussel_clam_integration_plan.md](mussel_clam_integration_plan.md).
It uses the same screen as the oyster backlog in
[candidate_studies.md](candidate_studies.md). Every accession, run count and
group structure below was read from ENA's deposited run and sample metadata on
**2026-10-12**, not from paper text. Publications were matched through PubMed on
the same date; a project with no publication listed had none that could be
found.

This is a working backlog, not a commitment. Claim a study with
`aree start-study` as for oysters ([parallel_study_runs.md](parallel_study_runs.md)),
and update this page when it is curated or rejected.

## How the screen was done

1. Every public RNA-Seq run in ENA for the six registered mussel and clam taxa
   was listed by project (`scripts/list_ena_rnaseq_projects.py`): **186
   projects**, from 8 for *M. trossulus* to 75 for Manila clam.
2. Every title was read. About 40 had a resilience-relevant design and enough
   runs to consider.
3. For those, the deposited design was recovered from sample attributes
   (`scripts/inspect_ena_designs.py`): groups and their sizes, the species
   recorded on each run, read depth, and layout.
4. Publications were matched for the strongest candidates, to confirm what
   phenotype was measured and whether two projects share animals.

## What decides the order

The oyster backlog's criteria carry over: raw reads public, at least three
biological replicates per group confirmed from run metadata, a measured
resilience outcome preferred over exposure only, and RNA-seq first. Two more
matter for these species:

- **Pooling needs the same species and the same phenotype.** Evidence pools on
  `species_taxid` + `phenotype` + `feature_type`. *M. edulis* and *M.
  galloprovincialis* never pool with each other, and an acidification study
  scored on larval viability never pools with one scored as acidification
  tolerance. A study that would be the only one in its species-and-phenotype
  group is worth less now than one that completes a pair.
- **Two projects from one lab are not automatically two studies.** Several
  candidates below come from the same group, sometimes from the same rearing
  program. Each pair is marked with what is known about whether the animals
  differ. Confirm from the methods before counting them as independent.

## Tier 1: the first within-species pools

The acceptance target for WP5 is at least one within-species pool with k ≥ 2
and real standard errors for mussels and for clams. Three pools are within
reach.

### A. Manila clam, heat (k = 2) — curate first

> **Claimed 2026-10-12** as `TUCCI2025_HEAT_PRIMING` (PRJNA1132140, draft PR
> #41), registered but **not runnable yet**. Reading the paper's methods found
> three things the deposit does not show:
>
> - **The libraries are 3′ Tag-Seq** (Lexogen QuantSeq, 75 bp single-end),
>   not full-length RNA-seq. The RNA-seq workflow's length correction (Salmon
>   effective length, tximport offsets in DESeq2) is wrong for 3′ counts, so
>   the workflow needs a 3′ Tag-seq mode first. `PRJNA906528` is from the same
>   group at the same depth and layout; check its library type before claiming.
> - **Tank is the unit of treatment** (two 25 L tanks per condition) and tank
>   membership is not deposited, so it cannot be fitted on tanks without the
>   authors' tank IDs.
> - **The heatwave clams are survivors.** RNA was taken 15 days after the
>   heatwave, from the digestive gland, after about 25% of naive clams had died.

| Project | Runs used | Design (from deposited metadata) | Publication | Size |
|---|---|---|---|---|
| `PRJNA1132140` | 40 RNA-seq, single-end, about 7 M reads | Primed vs naive × heatwave vs no heatwave, **n = 10 per group** | Tucci et al. 2025, *iScience* 28:113108, [10.1016/j.isci.2025.113108](https://doi.org/10.1016/j.isci.2025.113108) | 13.8 GB |
| `PRJNA906528` | 40 RNA-seq, single-end, about 7 M reads | 30 °C vs 20 °C for 30 days, digestive gland, **n = 20 per group** | none found | 7.9 GB |

Why first: the cheapest real pool in the screen (about 22 GB, small enough for
an 8 GB workstation), large n, and a measured outcome. Tucci et al. report
significantly higher survival in heat-primed clams under lethal heatwaves, so
primed vs naive under a heatwave is a contrast between groups that differ in a
measured survival outcome. That is the closest thing in this screen to the
oyster backlog's resistant-vs-susceptible lineages, though it is induced rather
than heritable, and survival is a group property, not measured per sequenced
animal. Register it as `thermal_tolerance` with the classification argued in
the YAML.

Curation notes:

- **Both projects deposit 16S amplicon runs from the same animals under the
  same attributes** (40 and 39 runs). Counting them doubles every group. Build
  the sample sheets with `--library-strategy RNA-Seq` (added for this case;
  the per-strategy breakdown is printed whenever strategies are mixed).
- **Pooling contrast.** The comparison that pools across the two is heatwave vs
  control in naive clams (`PRJNA1132140`) with 30 °C vs 20 °C (`PRJNA906528`).
  Make those the `meta_analysis_primary` comparisons; primed vs naive under a
  heatwave is a second, unpooled comparison.
- **Sex is unbalanced in `PRJNA906528`** (20 °C: 15 male, 5 female; 30 °C: 10
  and 10). Carry `sex` into the sheet and fit it as a covariate.
- **Independence.** Both come from the University of Padova group (the
  priming project lists Milan, Peruzza and Bargelloni; the 2023 project is
  unpublished). The 2023 project is a 30-day non-lethal exposure deposited
  before the priming experiment was published, so separate animals are likely,
  but confirm sampling dates or source cohort before counting k = 2.
- Tissue for `PRJNA1132140` is not in the deposit; the paper gives digestive gland.

```bash
aree fetch-samplesheet --bioproject PRJNA1132140 --study TUCCI2025_RPHIL_PRIMING \
  --condition-attribute Treatment --condition-attribute "Thermal stress" \
  --attribute Clam_ID --library-strategy RNA-Seq
aree fetch-samplesheet --bioproject PRJNA906528 --study <STUDY_ID> \
  --condition-attribute temp --attribute sex --library-strategy RNA-Seq
```

### B. Hard clam, salinity (k = 3)

| Project | Runs | Design | Publication | Size |
|---|---|---|---|---|
| `PRJNA975821` | 21 PE, 22 M reads | Gill. Acute and chronic hypo- (10 psu) and hypersalinity (45 psu) vs 30 psu, **n = 3 per group**; aliases AL1/AL5 and AH1/AH5 suggest two acute time points | Zhou et al. 2023, *Mar Environ Res* 192:106198, [10.1016/j.marenvres.2023.106198](https://doi.org/10.1016/j.marenvres.2023.106198) | 52.5 GB |
| `PRJNA1159784` | 9 PE | Gill. Three months at 15 or 35 ppt vs 25 ppt, **n = 3** | Zeng et al. 2024, *Fish Shellfish Immunol* 154:109888, [10.1016/j.fsi.2024.109888](https://doi.org/10.1016/j.fsi.2024.109888) | 5.2 GB |
| `PRJNA1288005` | 15 PE, 44 M reads | Hemocytes. 72 h at 5, 15, 25 (control), 35, 45 ppt, **n = 3** | Zeng et al. 2025, *Fish Shellfish Immunol* 166:110609, [10.1016/j.fsi.2025.110609](https://doi.org/10.1016/j.fsi.2025.110609) | 102 GB |

Three independent experiments from two labs (IOCAS Qingdao; University of
Florida), two of them on gill. All are `exposure_only` for transcription:
Zeng et al. 2025 measured survival (none died at 5 to 35 ppt; all died at
45 ppt by day 7) but per group, not per sequenced animal. Prespecify one
hyposalinity comparison per study as `meta_analysis_primary` (for example
chronic 10 psu vs 30 psu; 15 vs 25 ppt; 15 vs 25 ppt) so the pool compares
like with like. Only the needed runs of `PRJNA1288005` and `PRJNA975821` have to
be staged. `PRJNA1159784`'s run report carries no read counts; confirm depth
after download.

### C. *Mytilus edulis*, heat (k = 2, up to 4)

| Project | Runs | Design | Publication | Size |
|---|---|---|---|---|
| `PRJEB48959` | 35 PE, 51 M reads | Gill. 7 salinity × temperature groups, **n = 5** (ArrayExpress E-MTAB-11178). Control 23 ‰ at 5 °C; 30 and 33 °C at 23 ‰ and 15 ‰; 5 ‰ at 5 °C | Barrett et al. 2022, *Genes* 13:155, [10.3390/genes13010155](https://doi.org/10.3390/genes13010155) | 109 GB |
| `PRJEB41447` | 24 PE, 48 M reads | Outer-fjord animals warmed to 22 °C (n = 5) or 32 °C (n = 4) vs unwarmed (n = 5), plus field collections on a warm day (E-MTAB-9786) | Clark et al. 2021, *Sci Total Environ* 767:144366, [10.1016/j.scitotenv.2020.144366](https://doi.org/10.1016/j.scitotenv.2020.144366) | 85 GB |
| `PRJNA1050756` | 58 PE, 56 M reads | Two Atlantic Canada sources (Sober Island, NS; St. Peters Bay, PEI) × 20, 25, 30 °C; 30 °C has only 11 animals in total | none found | 319 GB |
| `PRJNA1377038` | 118 PE, 63 M reads | Diploid and triploid, slow and fast growers, 18 vs 26 °C | none found | 654 GB |

The two Greenland studies form the pool: same research group (British
Antarctic Survey and Aarhus), separate published experiments, comparable
contrasts (32 or 33 °C vs a cold control). Stage only the contrast runs (about
15 libraries from `PRJEB48959`, 9 from `PRJEB41447`).

Curation notes:

- **Species identity.** None of these deposits reports genotyping the animals,
  and *M. trossulus* co-occurs with *M. edulis* in Atlantic Canada. Every *M.
  edulis* study needs `species_identity` and, short of genotyping,
  `species_assignment_unverified` (validation enforces this).
- `PRJNA1050756`: 11 animals at 30 °C against 23 and 24 at the lower
  temperatures suggests mortality at 30 °C. If so, the 30 °C transcriptomes are
  survivors, which biases that arm. Check before using it.
- `PRJNA1377038` is the only *M. edulis* candidate with a growth phenotype
  (fast vs slow growers) under heat, which makes it the most resilience-relevant.
  At 654 GB it is a Tier 1 study to run on the larger workstation, not first.

## Tier 2: good, with a specific caveat

| Project | Species | Design | Context | Caveat |
|---|---|---|---|---|
| `PRJNA776346` | *M. galloprovincialis* | Larvae, pH 7.4 vs 8.1, static and at mid-trochophore, 3 replicate cultures | Acidification; larval abnormality measured (Kapsenberg et al. 2018, described in the deposit) | Cultures are pools of embryos from nine parental pairs, so n = 3 cultures. Alone in its *M. galloprovincialis* larval group, but the natural cross-species partner for the oyster larval OA studies once orthology exists (WP6). |
| `PRJNA954971` | *M. galloprovincialis* | 30 vs 18 °C, gill and mantle, n = 3 | Heat | Clean but small, and the only *M. galloprovincialis* heat study with n ≥ 3, so k = 1 for now. |
| `PRJNA611904` | *M. galloprovincialis* | Carbamazepine × low pH, n = 5, single-end | Acidification × pollutant | Low pH vs control is usable; the pharmaceutical arms are a different stressor. |
| `PRJNA688298` | *M. californianus* | Copper-exposed larvae sorted into normal vs abnormal, n = 3 pools per dose | Pollutant; within-exposure tolerant vs sensitive | Rare design (outcome measured on the sequenced animals), but a pollutant stressor and pooled larvae. |
| `PRJNA738278` | Manila clam | Two populations × −1, 12, 28 °C, n = 3 | Heat and cold, population contrast | Small n per cell; would add a third Manila clam heat study. |
| `PRJNA733332` | Manila clam | "Zebra 2" strain vs ordinary clams × ammonia, n = 3 | Ammonia stress, strain contrast | Selected-strain contrast is valuable; ammonia maps to `pollutant_exposure`. |
| `PRJNA621633` | Manila clam | Air exposure, re-immersion, control, n = 3 | `emersion_desiccation` | The first use of that stressor term; groups recoverable only from sample aliases. |
| `PRJNA1391265` | Manila clam | Chemical, thermal and mixed priming, 108 RNA-seq runs plus 108 amplicon | Heat priming | Design not in the sample attributes; recover it from the publication or submitter. 290 GB. |
| `PRJNA1289209` | Manila clam | Eggs from broodstock primed at low pH vs ambient, 30 runs | Acidification, parental priming | Larval outcomes measured per the description; design fields not deposited. |
| `PRJEB91723`, `PRJEB91724` | *M. edulis* | Noise exposure, then *Vibrio splendidus*, 24 runs each | Pathogen challenge | Infection arm n = 6 vs 18 controls across noise groups; two projects to reconcile. |
| `PRJNA1149726` | *M. edulis* | *Francisella halioticida* low- vs high-virulence strain vs sterile broth, hemocytes, 24 h | Pathogen challenge | Libraries are pools of animals, 4 per strain. |
| `PRJNA914271` | *M. edulis* | Hypoxia, low salinity, food × two sites × adult and seed, 112 runs | Multi-stressor | Design is in neither attributes nor aliases; 343 GB. |
| `PRJNA748538` | Manila clam | Low salinity vs control at 12 h, n = 3 | Salinity | Aliases (C123, T456) suggest each library pools three animals. |

## Same program, not independent: hard clam acidification

`PRJNA887747` (larvae to 10 months, acidified vs normal seawater, n = 4 per age),
`PRJNA900255` (extrapallial fluid and hemolymph after one year, n = 8 per tissue;
Schwaner et al. 2022, *Int J Mol Sci* 23:16066,
[10.3390/ijms232416066](https://doi.org/10.3390/ijms232416066)) and
`PRJNA909768` (mantle at 22 months, varying pCO2) all come from the Stony Brook
long-term acidification program. They look like three sampling points of one
rearing experiment, not three studies. Register them as **one** study with a
prespecified primary comparison, not as k = 3, unless the methods show separate
cohorts.

## Rejected, and why

| Project | Reason |
|---|---|
| `PRJNA824625` (*M. galloprovincialis*, salinity) | Libraries are pools per tank, two tanks per arm: the unit of replication is the tank, n = 2. |
| `PRJNA770496` (*M. galloprovincialis*, heat) | One library per temperature. |
| `PRJNA995578` (mussel and oyster mantle, acidification) | 12 runs per species deposited under one BioSample each: the runs are not shown to be distinct animals. |
| `PRJNA1285511`, `PRJNA933718` (Manila clam, heat) | 21 runs under one BioSample; 6 runs under two. Runs are not animals. |
| `PRJNA478917` (Manila clam, hypoxia) | One library per condition. |
| `PRJNA764366`, `PRJNA685920` (hard clam, heat × hypoxia; air exposure) | One library per condition or time point. |
| `PRJNA339676` (Manila clam, low salinity) | Two libraries per group. |
| `PRJNA1013758` (Manila and grooved carpet clam, heat × habitat) | 11 Manila clam libraries over four groups; the other species is not registered. |
| `PRJNA560413`, `PRJNA249058` | Population and phylogenomic sampling, no treatment. |
| `PRJNA1227433` (mussel, cockle, grooved carpet clam × *Vibrio*) | Three vs three per species at 11 M reads; kept in mind for WP6, not for pooling. |

## Not RNA-seq evidence, noted for later

- Bivalve transmissible neoplasia: `PRJNA749800` (*M. trossulus* BTN, 140 GB),
  `PRJNA1249915` (hard clam neoplasia). `disease_associated` only; see the plan.
- Byssal gland surveys (`PRJNA1225934`, *M. edulis*) and the multi-tissue atlases
  are annotation resources, not contrasts.
- *M. californianus* has the only sizeable mussel bisulfite deposits (94 runs,
  two projects); they are the first mussel methylation candidates once a mussel
  RNA-seq study is through.

## Sourcing notes

- ENA was queried by `tax_tree(<taxid>)` with `library_strategy="RNA-Seq"`.
  Mislabelled strategies are possible in both directions (see
  `LUTIER2022_OA_TIPPING`), so a project missing here may still hold RNA-seq.
- The species on every run was read from ENA, not assumed from the title:
  `PRJNA588016` (Gulf of Maine) mixes 20 *M. edulis* and 4 *M. trossulus* runs,
  and `PRJNA743680` (hard clam breeding) is half *Mercenaria campechiensis*.
- No ENA study record for these candidates carries a PubMed link; the
  publications above were found by searching PubMed for the design.
