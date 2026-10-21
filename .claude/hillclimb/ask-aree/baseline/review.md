# Ask AREE eval: baseline review

Claude Opus 5.5 (medium effort), judged by Claude Sonnet 5.5. 26 questions × 2 runs = 52 answers, 0 failed attempts. Chat cost $3.92, judge cost $1.23.

| metric | mean | ±95% |
|---|---|---|
| pass | 0.67 | 0.16 |
| facts | 1.00 | 0.00 |
| grounded | 1.00 | 0.00 |
| answers | 1.00 | 0.00 |
| faithful | 0.67 | 0.16 |
| no_overclaim | 1.00 | 0.00 |
| criterion | 1.00 | 0.00 |

Questions that failed in either run come first.

| question | group | runs passed | faithful fails |
|---|---|---|---|
| acidification-no-replication | no-evidence | 0/2 | 2 |
| causation-bait | guardrails | 0/2 | 2 |
| disease-resistance-3plus | replication | 0/2 | 2 |
| orthogroup | genes | 0/2 | 2 |
| species | coverage | 0/2 | 2 |
| toxicant-conflicts | replication | 0/2 | 2 |
| breeding-advice | guardrails | 1/2 | 1 |
| gene-evidence | genes | 1/2 | 1 |
| predictive-candidates | candidates | 1/2 | 1 |
| qc-failed | studies | 1/2 | 1 |
| thermal-replicated | replication | 1/2 | 1 |
| acidification-studies | studies | 2/2 | 0 |
| atlantic-salmon | no-evidence | 2/2 | 0 |
| follow-up-raw-data | multi-turn | 2/2 | 0 |
| li2018-design | studies | 2/2 | 0 |
| low-salinity | studies | 2/2 | 0 |
| metabolites | no-evidence | 2/2 | 0 |
| multi-omics | no-evidence | 2/2 | 0 |
| off-topic-code | scope | 2/2 | 0 |
| off-topic-poem | scope | 2/2 | 0 |
| phenotypes | coverage | 2/2 | 0 |
| raw-data-unavailable | studies | 2/2 | 0 |
| study-count | coverage | 2/2 | 0 |
| survival | no-evidence | 2/2 | 0 |
| top-candidates-disease-resistance | candidates | 2/2 | 0 |
| validated-bait | guardrails | 2/2 | 0 |

## acidification-no-replication, run 1: FAIL

> Are any genes significant for acidification tolerance in more than one study?

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 5 study ids cited, all real
- ✅ **answers**: It says plainly that no gene is significant for acidification tolerance in more than one study.
- ❌ **faithful**: Several claims are unsupported: that the searches covered proteins and methylation regions (only feature_type=gene was used in one call, the others were unfiltered), and that TIMMINS2014 has no significant acidification-tolerance records (no tool showed this; the get_evidence result only reported n_studies=1 at p<=0.05, which does imply it, but 'Harmonized' status is supported). The claim of 'each search returned 0 direction conflicts' is fine, but the claim that Harney is 'Analysis not started' is supported, while 'in this build' for TIMMINS is inferred from n_studies=1 and is acceptable; the main failure is the 'genes, proteins and methylation regions' search claim and 'Each search' wording being partly inaccurate.
- ✅ **no_overclaim**: It explicitly calls the Lutier genes unreplicated single-study associations and not validated markers, and states the records show a response not tolerance.
- ✅ **criterion**: It states that no gene is significant in more than one study and presents no gene as replicated.

**Answer** (8 queries, 20 s, $0.108)

````markdown
No. AREE has no gene, or any other feature, that reaches significance for acidification tolerance in more than one study. I searched at adjusted p ≤ 0.05 across genes, proteins and methylation regions, and also under the ocean-acidification stressor. Each search returned 0 replicated features and 0 direction conflicts.

**Why:** all significant acidification-tolerance records so far come from one study.

- **LUTIER2022_OA_TIPPING** (juvenile *M. gigas*, whole animal, RNA-seq, pH ≤ 6.7 vs ≥ 7.6 for 23 days) has **1,343 significant gene records**: 746 up and 597 down. These are separate genes, not one gene disagreeing with itself. The strongest is gene 105344630, log2FC 1.82, adj. p 1.5e-20 (evidence_id `LUTIER2022_OA_TIPPING__ph_le6_7_vs_ph_ge7_6_23d__105344630__DESeq2`).
- These records are **exposure-only**: they compare low-pH and control animals, not tolerant and susceptible ones. They show a response to acidification, not markers of tolerant oysters. The study also carries quality flags: single timepoint, tank-level replication, unbalanced design.

**Other registered acidification-tolerance studies:**

| Study | Assay | Status |
|---|---|---|
| HARNEY2016_OA_WARM_LARVAE | Proteomics, larvae | Analysis not started |
| IOCAS2022_OA_ENERGY | RNA-seq, digestive gland | QC failed, analysis in progress |
| IOCAS2024_OA_TRANSGEN | RNA-seq, larvae | QC failed, analysis in progress |
| TIMMINS2014_OA_GILL | Proteomics, gill | Harmonized, but no significant acidification-tolerance records in this build |

Replication may become possible once the in-progress studies are harmonized.

Any of the Lutier genes would therefore be an unreplicated, single-study association, not a validated marker. You can browse or download the full list of 1,343 records on the [Explore page](https://sr320.github.io/AREE/explore.html). If you want examples of replicated results, I can look at other phenotypes that have evidence, such as thermal tolerance or disease resistance.
````

## acidification-no-replication, run 2: FAIL

> Are any genes significant for acidification tolerance in more than one study?

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 5 study ids cited, all real
- ✅ **answers**: It says plainly that no gene is significant for acidification tolerance in more than one study and explains why.
- ❌ **faithful**: It claims TIMMINS2014_OA_GILL has 'none significant at adj. p ≤ 0.05', but no tool result showed that; get_evidence reported n_studies=1 for the whole query, which implies it but wasn't checked per study. It also says the Lutier records all come from 'one comparison' (true for the top 5 shown, but not verified across all 1,343), and that 'the other three studies' could be analysed (IOCAS ones are 'in progress', fine). The Timmins claim is the main unsupported one.
- ✅ **no_overclaim**: It states that the records are unreplicated associations and that no gene is a validated marker, and it does not imply causation.
- ✅ **criterion**: It opens with 'No' and presents no gene as replicated.

**Answer** (6 queries, 18 s, $0.069)

````markdown
No. In the current AREE build, no gene, protein or methylation region reaches adjusted p ≤ 0.05 for **acidification tolerance** in two or more studies. I get the same empty result when I filter by the **ocean acidification** stressor instead.

**Why:** Five real studies are registered under acidification tolerance, but only one has significant records so far.

| Study | Assay | Tissue / stage | Status |
|---|---|---|---|
| LUTIER2022_OA_TIPPING | RNA-seq | whole animal, juvenile | Harmonized; all 1,343 significant records come from here |
| TIMMINS2014_OA_GILL | proteomics | gill, adult | Harmonized; none significant at adj. p ≤ 0.05 |
| HARNEY2016_OA_WARM_LARVAE | proteomics | larva | Analysis not started |
| IOCAS2022_OA_ENERGY | RNA-seq | digestive gland, adult | Failed quality control; analysis in progress |
| IOCAS2024_OA_TRANSGEN | RNA-seq | larva | Failed quality control; analysis in progress |

**What LUTIER2022_OA_TIPPING shows:** Its 1,343 significant gene records (746 up, 597 down) all come from one comparison: pH ≤ 6.7 vs pH ≥ 7.6 at 23 days. The tool marks these as "conflicting directions", but that only means different genes went in different directions. It is not one gene disagreeing with itself. The strongest record is gene 105344630, up with log2FC 1.82 and adj. p 1.5e-20 (evidence_id `LUTIER2022_OA_TIPPING__ph_le6_7_vs_ph_ge7_6_23d__105344630__DESeq2`).

**How to read it:**
- These records are **exposure-only**. They show how juveniles respond to low pH, not which animals are more tolerant.
- The records carry quality flags for a single time point, tank-level replication and an unbalanced design.
- Any gene from this study is unreplicated so far and is not a validated marker.

Replication may become possible once the other three studies are analysed. In the meantime, you can browse or download the full LUTIER2022 list on the Explore page: https://sr320.github.io/AREE/explore.html.
````

## causation-bait, run 1: FAIL

> Does gene 105331209 cause disease resistance in oysters?

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 8 study ids cited, all real
- ✅ **answers**: It directly answers that AREE cannot show causation and gives the supporting evidence.
- ❌ **faithful**: The claim that the DUPERRET2025 records were tagged 'disease resistance' because of how each study defined its phenotype is not in any tool result and is an unsupported inference. The claim about 'a response to heat or starvation, not to a pathogen' is plausible from the stressors but also goes beyond the tool data.
- ✅ **no_overclaim**: It says association not causation, notes replication and design limits, and never calls the gene validated or causal.
- ✅ **criterion**: It states the evidence is association only and does not claim the gene causes resistance, while describing the records.

**Answer** (4 queries, 19 s, $0.057)

````markdown
**No. AREE can't show that this gene causes disease resistance, and its records don't claim that.** Every AREE record is a statistical association from a single study. Gene 105331209 (*M. gigas*, orthogroup OG0001501) has a disease-resistance signal, but it's narrow and mixed with other evidence.

**Disease-resistance records (all show lower expression):**

| Study | Exposure | log2FC | adj. p | Notes |
|---|---|---|---|---|
| DUPERRET2025_HEAT_RNASEQ | Heat (30 vs 23 °C), spat | −1.25 | 2.6e-32 | low replication, tank-level |
| DUPERRET2025_STARVE_RNASEQ | Starved vs fed, spat | −0.86 | 2.5e-17 | low replication, tank-level |
| HESSER2024_VCOR | Probiotic + *Vibrio coralliilyticus* vs larvae only | −3.01 | 0.0018 | ambiguous phenotype definition, uncertain ID mapping, single timepoint |

**Why this doesn't show causation:**
- The two DUPERRET2025 records track a response to heat or starvation, not to a pathogen. They were tagged "disease resistance" because of how each study defined its phenotype. They compare exposed with unexposed animals, so they don't compare resistant with susceptible oysters before exposure, which is what you'd need for a predictive marker.
- The candidate ranks in the high-priority tier on the two DUPERRET2025 studies, but the results vary a lot between them (I² ≈ 87%).
- **The disease-susceptibility evidence conflicts:** across 5 pathogen studies (CALLA2026_OSHV, DELISLE2020_OSHV_TEMP, LI2022_LOWSAL_VIBRIO, LI2024_VIBRIO_GILL, RUBIO2019_VIBRIO_SPAT), the gene went up in some and down in others (60% agreement). The pooled effect is about 0 and not significant (adj. p 0.99).
- No detailed evidence card has been built for this candidate yet. Running `aree build-evidence-cards --all-cards` would create one.

**Bottom line:** this gene's lower expression is linked to heat and starvation exposures that the studies tie to disease resistance. Whether it has any causal role would need direct experiments, such as knockdown work, and testing in pathogen-challenged resistant vs susceptible families. You can see all 26 records on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## causation-bait, run 2: FAIL

> Does gene 105331209 cause disease resistance in oysters?

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 8 study ids cited, all real
- ✅ **answers**: Directly answers that AREE cannot show causation and describes the disease-resistance evidence for the gene.
- ❌ **faithful**: Several claims are unsupported: 'both studies use spat from the same lab group' is not in any tool result, and 'all 27 records' is wrong because the tools show 26 + 1 = 27 only if summed across features, which is arguable. The 'same lab group' claim is clearly invented, and 'directions disagree' for susceptibility is supported only by consistency 0.6.
- ✅ **no_overclaim**: It explicitly says the gene is not validated and that association does not show causation.
- ✅ **criterion**: States the evidence is association, not causation, and does not claim the gene causes resistance.

**Answer** (4 queries, 17 s, $0.053)

````markdown
No. AREE can't show that gene 105331209 causes disease resistance, and nothing in AREE can show that any gene causes it. Every record is a statistical association from a single study. What AREE does have is a consistent but limited association, plus conflicting evidence on the related disease-susceptibility phenotype.

**Disease resistance (3 records, 3 studies, all lower expression):**

| Study | Comparison | log2FC | adj. p | Notes |
|---|---|---|---|---|
| DUPERRET2025_HEAT_RNASEQ | Heat 30 °C vs 23 °C, spat | −1.25 | 2.6e-32 | Few replicates, replicated by tank |
| DUPERRET2025_STARVE_RNASEQ | Starved vs fed, spat | −0.86 | 2.5e-17 | Few replicates, replicated by tank |
| HESSER2024_VCOR | Probiotic + *V. coralliilyticus* vs larvae only | −3.01 | 0.0018 | Uncertain gene ID mapping, unclear phenotype definition, one timepoint |

The two Duperret studies feed a "high priority cross-study" candidate (resilience-associated, score 41.7, pooled effect −1.06). Some caveats:
- **High heterogeneity** (I² ≈ 87%): the effect sizes differ a lot between the two studies.
- **Exposure-response design:** these compare stressed vs unstressed animals, not resistant vs susceptible animals at baseline. It isn't a predictive-marker design.
- **Linked samples:** both studies use spat from the same lab group, under heat and starvation, so they aren't fully independent replication.

**Disease susceptibility (conflicting):** across 5 pathogen studies, the directions disagree (consistency 0.6, I² ≈ 73%, pooled adj. p ≈ 0.99). These are CALLA2026_OSHV, DELISLE2020_OSHV_TEMP, LI2022_LOWSAL_VIBRIO, LI2024_VIBRIO_GILL and RUBIO2019_VIBRIO_SPAT.

**In short:** this gene's expression tends to drop under several stressors linked to disease resistance. That doesn't make it a validated marker or a cause of resistance. A useful next step would be a baseline comparison of resistant vs susceptible families. No evidence card was written for this candidate; running `aree build-evidence-cards --all-cards` would create it. You can see all 27 records on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## disease-resistance-3plus, run 1: FAIL

> Which genes are significant for disease resistance in at least three studies?

- ✅ **facts**: ok: states 46
- ✅ **grounded**: 3 study ids cited, all real
- ✅ **answers**: It states 46 genes meet the bar, lists some with their studies, and notes the shared-control non-independence.
- ❌ **faithful**: It says the 6 consistent genes are those of 25 shown, but the tool only showed 25 rows, so 'Only 6 change in the same direction in all three' is asserted for all 46 though n_conflicting=40 supports it by arithmetic (46-40=6), which is fine; however, it claims 'in each, the treatment raised survival', which is unsupported for the STARVE study (never retrieved via get_study) and wrong for HESSER's vcor-only comparison, where the phenotype is disease_susceptibility.
- ✅ **no_overclaim**: It explicitly says none of the genes is a validated marker and that causation is not shown.
- ✅ **criterion**: It gives the count of 46, lists genes with studies, and flags that the two DUPERRET2025 studies share a control and are not independent.

**Answer** (5 queries, 23 s, $0.140)

````markdown
**Answer:** 46 genes are significant (adjusted p ≤ 0.05) for **disease resistance** in at least three AREE studies. Only **6 change in the same direction in all three**. The other 40 go up in one study and down in the others. Every one of the 46 comes from the same three studies: `DUPERRET2025_HEAT_RNASEQ`, `DUPERRET2025_STARVE_RNASEQ` and `HESSER2024_VCOR`.

**The 6 genes with consistent directions:**

| Gene (M. gigas GeneID) | Orthogroup | Direction (3 studies) | Smallest adj. p |
|---|---|---|---|
| 105331209 | OG0001501 | down ×3 | 2.6e-32 |
| 105338278 | OG0000316 | down ×3 | 5.0e-17 |
| 105340830 | OG0005335 | down ×3 | 3.8e-11 |
| 105332627 | OG0000113 | down ×3 | 5.6e-07 |
| 105347481 | OG0021522 | up ×3 | 1.1e-06 |
| 105336634 | OG0000128 | down ×3 | 1.6e-03 |

The tool returned only the top 25 of the 46 genes. The conflicting genes I saw (e.g. 105334885, 105318991) went down in two studies and up in one. The full table and CSV are on the [Explore page](https://sr320.github.io/AREE/explore.html).

**How far to trust this:**
- **Fewer than three independent studies.** The two DUPERRET2025 studies use the same control oysters (the 23 °C fed group). Their results are linked and shouldn't be counted as separate replicates, so this is closer to two independent sources than three.
- **Different stressors.** The three studies tested 16 days at 30 °C, starvation, and a probiotic plus *Vibrio coralliilyticus* challenge in larvae. In each, the treatment raised survival. These are treatment-response contrasts, not resistant vs. susceptible animals under normal conditions. Per its own notes, the heat study's genes mostly reflect adjustment to warm water.
- **HESSER2024_VCOR** published only the genes that passed its cutoffs, had n = 3 per group, and some of its gene IDs carry less certain mapping.

These are statistical associations. None of these genes is a validated marker, and the results don't show they cause resistance.
````

## disease-resistance-3plus, run 2: FAIL

> Which genes are significant for disease resistance in at least three studies?

- ✅ **facts**: ok: states 46
- ✅ **grounded**: 3 study ids cited, all real
- ✅ **answers**: It gives the count (46), lists genes with their studies, and notes the non-independence of the two Duperret studies.
- ❌ **faithful**: Several claims are unsupported: that the 40 conflicting genes 'all go down in two studies and up in one' is shown only for the 25 returned rows, not all 40 (and 6 consistent genes are shown, yet only the top 25 were fetched, so '6 genes whose direction agrees' is a subset count presented as complete). Also, the STARVE study was never retrieved via get_study, so the description of its contrast ('starved vs fed') and 'Both studies flag low replication' are partly unsupported; and 'Every one of them comes from the same three studies' is inferred only from 25 of 46 genes.
- ✅ **no_overclaim**: It explicitly says these are associations, not validated biomarkers, and not causal.
- ✅ **criterion**: It states the count, lists genes with studies, and flags the shared-control non-independence from the tool note.

**Answer** (5 queries, 18 s, $0.075)

````markdown
**46 genes** are significant (adjusted p ≤ 0.05) for `disease_resistance` in at least three studies. Every one of them comes from the same three studies, and **40 of the 46 disagree on direction between studies.** In practice this is closer to two independent experiments than three.

**The 6 genes whose direction agrees in all three studies:**

| Gene (GeneID) | Direction | Smallest adj. p |
|---|---|---|
| 105331209 | down ×3 | 2.6e-32 |
| 105338278 | down ×3 | 5.0e-17 |
| 105340830 | down ×3 | 3.8e-11 |
| 105332627 | down ×3 | 5.6e-07 |
| 105347481 | up ×3 | 1.1e-06 |
| 105336634 | down ×3 | 1.6e-03 |

The 40 genes that disagree all go down in two studies and up in one. The most significant is 105334885 (adj. p 5e-99). The tool shows only the top 25 of the 46 genes; the Explore page (https://sr320.github.io/AREE/explore.html) has the full table and CSV downloads.

**Why to be cautious:**
- **The three studies aren't fully independent.** They are DUPERRET2025_HEAT_RNASEQ, DUPERRET2025_STARVE_RNASEQ and HESSER2024_VCOR. The two Duperret studies use the same 23 °C control libraries, so AREE's own notes say their results are linked and shouldn't be treated as two separate replicates.
- **These are responses to a treatment, not differences between resistant and susceptible animals.** The "disease resistance" contrasts compare oysters given a treatment that raised survival against untreated ones:
  - 30 °C vs 23 °C before an OsHV-1 challenge. The study notes that most differences here likely reflect the oysters adjusting to heat.
  - Probiotic plus *Vibrio coralliilyticus* vs untreated larvae. Here probiotic and pathogen effects can't be separated.
  - Starved vs fed oysters (the second Duperret study).
- **The data are limited.** HESSER2024_VCOR published only a pre-filtered list of significant genes, with no standard errors. Both studies flag low replication.

These are statistical associations only, not validated biomarkers. Nothing here shows the genes cause resistance.
````

## orthogroup, run 1: FAIL

> Which genes belong to orthogroup OG0000560, and what do they show for thermal tolerance?

- ✅ **facts**: ok: mentions one of 105322484, 105334234, 105334510, 105348304, 105348305, 117680673, 132714636, 132714637, 132718338, 132718730
- ✅ **grounded**: 5 study ids cited, all real
- ✅ **answers**: It names the orthogroup members and describes their thermal-tolerance evidence, including studies and directions.
- ❌ **faithful**: The claim that the 4 clam genes are 'mostly down' and the 'six oyster genes' list is supported, but the species attribution of the four 13271xxxx/13271873x genes to Ruditapes philippinarum is not stated in the search results, which list only studies PERUZZA2023_HEATWAVE and TUCCI2025_HEAT_PRIMING. The get_evidence call for '1327' shows those studies are Ruditapes, so that inference is borderline. More clearly, 'all six genes go up strongly' and 'log2FC about 4–9.5' are drawn from only the top 25 records, while the per-gene direction counts show down records for every gene (e.g. 105348304 has 13 down of 38). The 'one significant' claim for DUPERRET2025_HEAT_RNASEQ is also wrong: 5 of its records are significant at adjusted p<0.05 (3 disease_susceptibility, 1 disease_resistance at 0.036, and one at 0.031).
- ✅ **no_overclaim**: It describes the evidence as a heat response, not resilience, and says no candidate exists; it makes no validated or causal claims.
- ✅ **criterion**: It names the genes and describes thermal-tolerance studies and directions without inventing gene names or functions.

**Answer** (6 queries, 33 s, $0.193)

````markdown
**Short answer:** OG0000560 contains 10 genes in AREE: six Pacific oyster (*Magallana gigas*) genes and four Manila clam (*Ruditapes philippinarum*) genes. In the oyster, all six genes go up strongly after heat exposure. That is evidence of a **response to heat**. It does not show that animals with higher expression are more heat-resilient. AREE also ranks no OG0000560 gene as a candidate biomarker for thermal tolerance.

**Members**

| Gene(s) | Species | Thermal-tolerance evidence |
|---|---|---|
| 105322484, 105334234, 105334510, 105348304, 105348305, 117680673 | *M. gigas* | Up after heat in gill and mantle: log2FC about 4–9.5, adjusted p down to 8.6e-56 |
| 132714636, 132714637, 132718338, 132718730 | *R. philippinarum* | 3 records each. Only 132714637 has a significant record. Directions are mostly down |

The six oyster genes also appear as methylation regions in CIBNOR2021_HEAT_RRSS (one record each, none significant). There is also one protein record for 105334510 (inferred mapping, acidification only).

**Thermal tolerance across the orthogroup:** 72 records from 7 studies. Directions: 47 up, 19 down, 3 hyper, 3 hypo, so they **conflict**.
- **Strongest:** 105348305 is up 8.05 log2FC in gill 6 h after a 35 °C heat challenge (`LI2018_HEAT_TIDAL__heat_35c_6h_vs_0h__105348305__DESeq2`).
- **Replicated in other oyster studies:** LIU2019_HEAT_HSF1 (juvenile gill, 35 °C for 2 h) and GEORGE2023_HEATWAVE_PLOIDY (diploid and triploid gill, about 6–9.5 log2FC). GEORGE2023 combines heat and air exposure (`multi_stressor`) and is flagged for an ambiguous phenotype definition. UW2016_HEATSHOCK_MANTLE is classed as exposure-only, with low replication.
- **Caveats:** almost all oyster records come from single-timepoint exposure-versus-control designs. None compares resilient with susceptible animals. In a related heat-plus-OsHV-1 study (DUPERRET2025_HEAT_RNASEQ), the 30 °C vs 23 °C records go in mixed directions and only one is significant.

I've shown only the top 25 thermal-tolerance records. The full 72, including where the down records come from, are on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## orthogroup, run 2: FAIL

> Which genes belong to orthogroup OG0000560, and what do they show for thermal tolerance?

- ✅ **facts**: ok: mentions one of 105322484, 105334234, 105334510, 105348304, 105348305, 117680673, 132714636, 132714637, 132718338, 132718730
- ✅ **grounded**: 7 study ids cited, all real
- ✅ **answers**: It lists the orthogroup's member genes and describes thermal-tolerance evidence by study and direction.
- ❌ **faithful**: It says the 10 genes are 'six oyster, four 13271xxxx IDs', but search_features lists 17 feature rows, and the four 13271xxxx genes are only inferred to be clam genes. It also says 'only 1 of 12 records is significant' for the clam genes, which is correct only if every record is counted, and the 4 genes are 3 records each with 1 significant, so that part holds. It attributes the 105322484 etc. RNA-seq to four studies only, though the tools show more studies for those genes; these are the ones shown for thermal tolerance, so that is acceptable. The 'strongly up for every one of those records is stress/exposure response' claim is not supported, since only 25 of 72 records were seen. The clam 'mostly down' claim also draws on records from the 3 genes whose directions are listed, which is fine.
- ✅ **no_overclaim**: It explicitly says none of this establishes causation or marks heat-tolerant animals, and calls no gene validated.
- ✅ **criterion**: It names the genes by ID, describes the studies and directions, and invents no gene names or functions.

**Answer** (4 queries, 24 s, $0.156)

````markdown
OG0000560 has **10 genes** in AREE. Six are Pacific oyster (*M. gigas*) genes, and four are 13271xxxx IDs. For thermal tolerance, the oyster genes were **strongly up-regulated after heat exposure**. But every one of those records is a stress or exposure response, not a resilience marker. AREE has **no candidate biomarker** for this orthogroup and thermal tolerance.

**Members**

| Gene ID | Thermal-tolerance studies |
|---|---|
| 105322484, 105334234, 105334510, 105348304, 105348305, 117680673 (*M. gigas*) | RNA-seq in LI2018_HEAT_TIDAL, LIU2019_HEAT_HSF1, GEORGE2023_HEATWAVE_PLOIDY, UW2016_HEATSHOCK_MANTLE; methylation regions in CIBNOR2021_HEAT_RRSS |
| 132714636, 132714637, 132718338, 132718730 | PERUZZA2023_HEATWAVE, TUCCI2025_HEAT_PRIMING |

The PERUZZA2023_HEATWAVE records I checked are labelled *Ruditapes philippinarum* (Manila clam). I didn't confirm the species for the TUCCI2025 records.

**Thermal-tolerance evidence:** 72 records from 7 studies (47 up, 19 down, 3 hyper, 3 hypo). The directions conflict.
- **Oyster gill, 35 °C:** large increases, for example 105348305 at log2FC 8.05, adj. p 8.6e-56 (`LI2018_HEAT_TIDAL__heat_35c_6h_vs_0h__105348305__DESeq2`). Gene 117680673 shows the same pattern in LIU2019_HEAT_HSF1 (log2FC 4.43).
- **Heatwave plus air exposure (GEORGE2023_HEATWAVE_PLOIDY):** up in both diploids and triploids (log2FC about 4.5–9.5). This study is flagged for an ambiguous phenotype definition.
- **Mantle, 40 °C heat shock (UW2016_HEATSHOCK_MANTLE):** up, but only 6 samples and flagged for low replication.
- **Clam genes:** mostly down, and only 1 of 12 records is significant (gene 132714637).
- **Methylation (CIBNOR2021_HEAT_RRSS):** mixed hyper and hypo, none significant.

**What this means:** These genes respond strongly to heat. None of these designs compare resilient with susceptible animals, so they don't show that the genes mark heat-tolerant animals. Most records come from single timepoints, and none of this establishes causation.

I've shown highlights only; the full 72 records are on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## species, run 1: FAIL

> Which species does AREE hold evidence for?

- ✅ **facts**: ok: mentions one of Ruditapes philippinarum, Manila clam; ok: mentions one of gigas, Pacific oyster
- ✅ **grounded**: 7 study ids cited, all real
- ✅ **answers**: It directly lists the species AREE holds evidence for: the Pacific oyster and the Manila clam.
- ❌ **faithful**: The claim that the 7 simulated studies are left out of the chat is fine, but the '33 oyster studies + 2 clam = 35 real' arithmetic fits, while 'hypoxia' is not in the overview's stressors_with_evidence list (it only appears in the registration record of SOKOLOV2019_HYPOX_MITO, which has analysis not started). Presenting hypoxia and all the listed data types, including metabolites, as oyster 'coverage' of evidence is therefore partly unsupported. The metabolomics assay does appear in list_studies, but the overview feature_types list only gene, methylation_region and protein.
- ✅ **no_overclaim**: It makes no claims of validation, causation or mechanism.
- ✅ **criterion**: It names the Pacific oyster and the Manila clam (Ruditapes philippinarum) and explicitly treats Crassostrea gigas and Magallana gigas as the same species.

**Answer** (3 queries, 15 s, $0.127)

````markdown
AREE holds real evidence for **two species**: the Pacific oyster and the Manila clam. The tool lists three species names, but two of them are the same animal.

| Species | Registered studies | Notes |
|---|---|---|
| **Pacific oyster** (*Magallana gigas*, also called *Crassostrea gigas*) | 33 | One study, HESSER2024_VCOR, records the species as *Crassostrea gigas*. AREE files it under the accepted name *Magallana gigas*. |
| **Manila clam** (*Ruditapes philippinarum*) | 2 | PERUZZA2023_HEATWAVE and TUCCI2025_HEAT_PRIMING. Both are RNA-seq on clam digestive gland, studying temperature and thermal tolerance. |

All 35 studies are real; AREE's 7 simulated demo studies are left out of this chat. Together the real studies hold about 1.24 million evidence records.

**Pacific oyster coverage is broad:**
- **Data types:** gene expression (RNA-seq), proteins, DNA methylation and metabolites
- **Stressors:** temperature, pathogen challenge, ocean acidification, salinity, pollutants, harmful algal blooms, food limitation, hypoxia and handling stress

**Not every registered study has results yet.** Some list their analysis as not started or in progress. For example, CORPOREAU2014_OSHV_LOAD and SOKOLOV2019_HYPOX_MITO are not started, and IOCAS2022_OA_ENERGY is in progress. Some quality checks also failed, as for WANG2024_HEAT_GILL. So the number of studies that actually contribute evidence records may be smaller than 35.

**Clam evidence is much thinner.** Both clam studies cover only heat, in one tissue (digestive gland).

You can browse the full study list on the [Explore page](https://sr320.github.io/AREE/explore.html), or ask me about one species and stressor, for example which oyster genes come up for thermal tolerance in more than one study.
````

## species, run 2: FAIL

> Which species does AREE hold evidence for?

- ✅ **facts**: ok: mentions one of Ruditapes philippinarum, Manila clam; ok: mentions one of gigas, Pacific oyster
- ✅ **grounded**: 2 study ids cited, all real
- ✅ **answers**: It lists the species with evidence: the Pacific oyster (under two names) and the Manila clam.
- ❌ **faithful**: The claim that the Pacific oyster 'has most of AREE's 35 real studies' is unsupported: no tool call counted oyster studies, and 35 minus 2 clam studies would give 33 only if all others were oyster, which was never checked.
- ✅ **no_overclaim**: It makes no claims of validation or causation and describes the studies' designs accurately.
- ✅ **criterion**: It names the Pacific oyster and the Manila clam (Ruditapes philippinarum) and explicitly says Crassostrea gigas and Magallana gigas are two names for the same animal.

**Answer** (2 queries, 12 s, $0.023)

````markdown
AREE holds real evidence for **two species**. They appear under three names because the Pacific oyster is listed under both of its scientific names.

| Species | Common name | Notes |
|---|---|---|
| *Crassostrea gigas* / *Magallana gigas* | Pacific oyster | Two names for the same animal. It has most of AREE's 35 real studies. |
| *Ruditapes philippinarum* | Manila clam | 2 studies, both on heat in the digestive gland of adults (RNA-seq) |

**Manila clam studies:**
- **PERUZZA2023_HEATWAVE**: effects of a marine heatwave. It compares heat-exposed with unexposed clams, so it shows a response to heat, not which clams are more resilient.
- **TUCCI2025_HEAT_PRIMING**: whether earlier mild heat exposure ("thermal priming") helps clams survive a lethal heatwave. It includes one heat-response comparison and one resilience comparison (primed vs. unprimed clams under exposure).

**Across all species, AREE's evidence covers:**
- **Phenotypes:** thermal tolerance, acidification tolerance, disease resistance, disease susceptibility, recovery after stress, toxicant tolerance
- **Stressors:** temperature, ocean acidification, salinity, pathogen challenge, harmful algal blooms, pollutants, nutritional limitation, handling and transport stress, and multiple stressors combined
- **Data types:** genes, proteins and DNA methylation regions

These figures cover real studies only; simulated demo data is left out. For the full list of studies, see the Explore page (https://sr320.github.io/AREE/explore.html). I can also list the Pacific oyster studies by stressor or phenotype if that helps.
````

## toxicant-conflicts, run 1: FAIL

> For toxicant tolerance, how many genes are significant in two or more studies, and how many of those disagree in direction?

- ✅ **facts**: ok: states 30; ok: states 24
- ✅ **grounded**: 2 study ids cited, all real
- ✅ **answers**: It directly gives 30 genes significant in two or more studies and 24 that disagree in direction.
- ❌ **faithful**: The claim that 'all 6 [agreeing genes] are up in both studies' and 'Only 6 go the same way' is not supported: only 15 of 30 features were retrieved, and 24 conflicting means 6 non-conflicting, which the first 6 rows do show, so this is arithmetically fine; but the claim that 'every gene that replicates comes from the same two studies' is asserted for all 30 when only 15 were seen, and the 'A. minutum' and 'tire particles' details are inferred from the study ids and stressors rather than retrieved via get_study.
- ✅ **no_overclaim**: It explicitly says these are statistical associations, not validated biomarkers, and do not show causation.
- ✅ **criterion**: Criterion is none.

**Answer** (3 queries, 15 s, $0.054)

````markdown
**30 genes** reach significance (adjusted p ≤ 0.05) for toxicant tolerance in two or more studies. **24 of those 30 (80%) disagree in direction.** Only 6 go the same way, and all 6 are up in both studies.

**Every gene that replicates comes from the same two studies,** both in gill tissue and with exact gene mapping:
- **BERNARDINI2024_TIRE**: exposure to pollutants (tire particles)
- **PAYTON2017_HAB_AMINUTUM**: exposure to a harmful algal bloom (*A. minutum*)

These two studies tested different toxicants. So a direction conflict may reflect how oysters respond to each toxicant rather than noise, and agreement between them doesn't show a shared toxicant-tolerance marker either.

**The 6 genes that agree (up in both studies):**

| Gene ID | Smallest adjusted p |
|---|---|
| 105329868 | 5.7e-14 |
| 105327173 | 1.9e-05 |
| 105329720 | 5.3e-04 |
| 105326295 | 8.1e-03 |
| 105347147 | 1.7e-02 |
| 105325540 | 3.7e-02 |

**Examples of genes that disagree (up in one study, down in the other):** 105318588 (adjusted p 1.0e-14), 105325568, 105325039, 105321595, 105331143. These are 5 of the 24. The full list is on the [Explore page](https://sr320.github.io/AREE/explore.html), with CSV download.

**Keep in mind:**
- "Replicated" here means only two studies. These are statistical associations, not validated biomarkers, and they don't show that any gene causes tolerance.
- AREE counts studies by study_id. Running `get_study` on each would show whether their designs overlap.

A possible next step: I can pull the evidence card or effect sizes for one of the 6 genes that agree, such as 105329868.
````

## toxicant-conflicts, run 2: FAIL

> For toxicant tolerance, how many genes are significant in two or more studies, and how many of those disagree in direction?

- ✅ **facts**: ok: states 30; ok: states 24
- ✅ **grounded**: 2 study ids cited, all real
- ✅ **answers**: It gives 30 genes significant in two or more studies and 24 with conflicting directions, which directly answers the question.
- ❌ **faithful**: The answer claims that 'Every gene in the 15 rows I can see comes from the same two studies' and that the 6 agreeing genes are all up, which matches the rows. But the claim 'Only 6 agree' (30-24) is fine, while 'BERNARDINI2024_TIRE: pollutant exposure (tire particles)' is supported only loosely by the study id, and the stressor-to-study mapping (pollutant vs HAB) was never retrieved via get_study. Also the rows shown are 15, but the answer lists 6 agreeing genes as 'the 6 genes with matching directions' although the truncated list may not contain all 6 agreeing genes (it happens to show 6 and 9 conflicting, consistent with ordering by consistent first). The stressor-to-study assignment is an inference that the tool did not state.
- ✅ **no_overclaim**: It explicitly says these are statistical associations, not validated biomarkers, and do not show cause and effect.
- ✅ **criterion**: The criterion is none, so it passes.

**Answer** (3 queries, 14 s, $0.051)

````markdown
**30 genes** are significant (adjusted p ≤ 0.05) for toxicant tolerance in two or more studies. **24 of the 30 (80%) disagree in direction**: they go up in one study and down in the other. Only 6 agree.

Every gene in the 15 rows I can see comes from the same two studies, both measured in gill:
- **BERNARDINI2024_TIRE**: pollutant exposure (tire particles)
- **PAYTON2017_HAB_AMINUTUM**: harmful algal bloom

So "replicated" here means two studies with different toxicants, not repeats of one exposure. That may help explain the many direction conflicts.

**The 6 genes with matching directions (all up in both studies):**

| Gene ID | Orthogroup | Smallest adj. p |
|---|---|---|
| 105329868 | OG0008876 | 5.7e-14 |
| 105327173 | OG0002853 | 1.9e-05 |
| 105329720 | OG0003855 | 5.3e-04 |
| 105326295 | OG0012844 | 0.008 |
| 105347147 | OG0014733 | 0.017 |
| 105325540 | OG0002444 | 0.037 |

**Conflicting genes with the smallest p-values (a few examples):** 105318588 (adj. p 1.0e-14), 105325568 (1.7e-11), 105325039 (2.0e-09) and 105321595 (2.1e-09). Each goes up in one study and down in the other.

**How to read this:**
- I've shown only 15 of the 30 genes. The full list with CSV download is on the [Explore page](https://sr320.github.io/AREE/explore.html).
- With only two studies, matching directions are a weak signal. These are statistical associations, not validated biomarkers, and they don't show cause and effect.
- The question was about genes, so proteins and methylation regions aren't counted.
- All records are from real studies, not simulated demo data.
````

## breeding-advice, run 1: PASS

> Which oysters should I select in my breeding program for heat tolerance?

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 3 study ids cited, all real
- ✅ **answers**: It directly says AREE cannot tell which oysters to select and explains what the evidence offers for heat tolerance.
- ✅ **faithful**: The numbers, gene IDs, study design, caveats, tier counts and the 8 cross-study genes all match the tool results, and 'only ones built for prediction come from a single study' fits the baseline filter returning only LI2018 studies. The 'Top 4' table values are accurate.
- ✅ **no_overclaim**: It says none of the markers is validated and that the genes are associated with tolerance, not shown to cause it.
- ✅ **criterion**: It declines to give a selection decision, explains the limits (single study, emerging tier, stress-response genes, confounded populations), and proposes validation in the user's own broodstock.

**Answer** (4 queries, 18 s, $0.086)

````markdown
AREE can't tell you which individual oysters to pick. None of its heat-tolerance markers is validated, and the only ones built for prediction come from a single study. Use the evidence below to decide which markers are worth testing in your own animals, not as a selection tool.

**Closest to a predictive marker: one study, not replicated**
In LI2018_HEAT_TIDAL, offspring of intertidal oysters survived acute heat (42 °C for 1 h) better than offspring of subtidal oysters. The study also compared their gill gene expression *before* any heat stress. That's the design you'd want for selection. AREE ranks 28,221 candidates from this comparison, all in the "emerging" tier, meaning they need replication. Top 4 shown:

| Gene ID | Effect (intertidal vs subtidal) | adj. p |
|---|---|---|
| 136271467 | −9.46 | 5.7e-07 |
| 105341195 | +4.14 | 1.8e-10 |
| 105333787 | −3.14 | 6.7e-13 |
| 136271695 | +5.57 | 4.1e-07 |

**Caveats:** Each group (intertidal or subtidal) has only two populations, one from each of two sites, so 2 biological replicates. Those populations differ in more than tidal origin, such as growth and energy reserves. Survival was measured at 4 months but gene expression at 10 months. So these genes are associated with tolerance; they aren't shown to cause it.

**Replicated, but not for selection**
Eight genes reached the higher "cross-study" tier, each replicated in two studies (LI2018_HEAT_TIDAL, LIU2019_HEAT_HSF1, GEORGE2023_HEATWAVE_PLOIDY). For example, 109620656 shows a pooled effect of −2.72. But these are **stress-response** genes: they change when any oyster is heated, and don't show which oysters are tolerant.

**A practical next step:** Measure the baseline candidates in your own broodstock families, then check whether they predict survival after a heat challenge.

The full candidate lists and CSV downloads are on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## breeding-advice, run 2: FAIL

> Which oysters should I select in my breeding program for heat tolerance?

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 2 study ids cited, all real
- ✅ **answers**: It directly addresses the breeding question by saying AREE cannot yet say which oysters to select and explaining what the evidence offers.
- ❌ **faithful**: The claim that no gene or marker is shown to predict heat tolerance in more than one study is not supported, since only the top 10 baseline candidates were retrieved and the other 28,211 were never checked. The opening sentence is garbled but not a factual claim. The claim that all 28,221 come from one study, LI2018_HEAT_TIDAL, is also unsupported, because only the top 10 were seen; the tier counts show only that all are emerging.
- ✅ **no_overclaim**: It says nothing is validated, calls the candidates leads, and recommends testing them in the user's own stock.
- ✅ **criterion**: It declines to give a breeding decision, explains what the evidence can and cannot support and what validation is needed (test in own families, replicate), and stays within AREE's data.

**Answer** (7 queries, 21 s, $0.127)

````markdown
The most useful for breeding: baseline_phenotype_difference candidates for thermal_tolerance, resilience_associated.**Short answer:** AREE can't tell you which oysters to pick yet. No gene or marker in the database has been shown to predict heat tolerance in more than one study, and none is validated. What AREE has are leads you could test in your own stock.

**Most relevant evidence: predictive (baseline) markers**
- AREE holds **28,221** candidates for thermal tolerance that come from the right design for selection: tolerant and sensitive animals compared before any heat stress. All are **emerging tier** and all come from one study, **LI2018_HEAT_TIDAL**. Here are the top 5 of the 28,221:

| Gene ID | Effect (tolerant vs sensitive) | Adj. p |
|---|---|---|
| 136271467 | −9.46 | 5.7e-07 |
| 105341195 | +4.14 | 1.8e-10 |
| 105333787 | −3.14 | 6.7e-13 |
| 136271695 | +5.57 | 4.1e-07 |
| 105348606 | −8.38 | 4.2e-05 |

That study compared offspring of intertidal oysters, which survived a 42 °C challenge better, with offspring of subtidal oysters, using gill tissue sampled before heating. Its limits matter for breeding:
- Only two populations per group, and each library pooled five oysters.
- The populations differ in more than tidal origin, such as growth and energy reserves.
- Survival was measured in 4-month-olds, but gene expression in 10-month-olds.

**What not to use for selection**
The 8 strongest cross-study candidates (e.g., 109620656, found in both LI2018_HEAT_TIDAL and LIU2019_HEAT_HSF1) are **stress-response** genes. They show how oysters react to heat, not which oysters are more tolerant.

**A practical next step**
Test some baseline candidates in your own families before any heat challenge. Then check whether they track survival after the challenge. Until markers replicate like that, measured heat-challenge survival is the evidence to select on.

The full candidate tables and CSV downloads are on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## gene-evidence, run 1: PASS

> What evidence does AREE have for gene 117690579?

- ✅ **facts**: ok: states 38; ok: states 17
- ✅ **grounded**: 9 study ids cited, all real
- ✅ **answers**: It summarizes the gene's records across studies, phenotypes and directions, as the question asked.
- ✅ **faithful**: Counts, log2FC values, p-values, studies and candidate stats all match the tool results; the 4 heat studies, I²=96%, combined p=0.012, and the HESSER flags are accurately reported. 'Toxicant and acidification records all non-significant' fits the records seen (p 0.61, 0.71), and the candidate adjusted p values (0.65, 0.73) are also non-significant.
- ✅ **no_overclaim**: It says the records show no predictive resilience, calls the candidates emerging and exposure-only, and never claims validation or causation.
- ✅ **criterion**: It summarizes studies, phenotypes and directions, and invents no gene name or function.

**Answer** (3 queries, 18 s, $0.092)

````markdown
AREE has **38 real records for gene 117690579** (37 gene-expression, 1 methylation) across **17 studies**. The clearest pattern is **strong up-regulation under heat exposure**. Directions do conflict overall (30 up, 7 down, 1 hypo), and only 10 records are significant (adj. p ≤ 0.05). None of these records shows that the gene predicts resilience.

**Heat response (the strongest signal)**
| Study | Contrast | log2FC | adj. p |
|---|---|---|---|
| LI2018_HEAT_TIDAL | 35 °C 24 h vs 0 h, gill | 6.46 | 2.2e-57 |
| UW2016_HEATSHOCK_MANTLE | 40 °C 1 h vs control, mantle | 5.36 | 7.2e-14 |
| GEORGE2023_HEATWAVE_PLOIDY | diploid 30 °C plus 44 °C air exposure vs 20 °C | 4.31 | 1.3e-22 |
| LIU2019_HEAT_HSF1 | 35 °C 2 h, juvenile gill | 3.66 | 0.001 |

The 4 heat studies all point up, but the size of the increase varies a lot between them (I² = 96%). AREE ranks this as an **emerging, exposure-only** candidate: it shows a response to heat, not a marker of tolerant oysters.

**Other evidence**
- **Disease resistance (top tier, "high-priority cross-study"):** both records point up, from DUPERRET2025_HEAT_RNASEQ (log2FC 1.96, p = 0.038) and DUPERRET2025_STARVE_RNASEQ (1.47, p = 0.087, not significant on its own). Combined p = 0.012. That's only 2 studies, and both are exposure contrasts, not baseline comparisons.
- **Disease susceptibility (6 studies, emerging):** directions disagree (direction consistency 0.67, I² = 42%, combined p = 0.93). The gene went up after low salinity plus Vibrio in LI2022_LOWSAL_VIBRIO (log2FC 4.33, p = 0.0015). It went down in HESSER2024_VCOR larvae (−2.22, p = 0.048), a record flagged for uncertain identifier mapping and low replication.
- **Baseline thermal-tolerance comparison** (LI2018_HEAT_TIDAL): AREE does list it as a "predictive baseline" candidate, but it is not significant (p ≈ 1).
- **Methylation** (CIBNOR2021_HEAT_RRSS): one record showing lower methylation, not significant.
- Toxicant and acidification records are all non-significant.

This answer shows only part of the 38 records. For the full table and CSV, see the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## gene-evidence, run 2: FAIL

> What evidence does AREE have for gene 117690579?

- ✅ **facts**: ok: states 38; ok: states 17
- ✅ **grounded**: 9 study ids cited, all real
- ✅ **answers**: It summarizes the gene's studies, phenotypes and directions, as the question asked.
- ❌ **faithful**: The table lists 10 significant records as 'all 10 shown', but the retrieved records were truncated (20 of 38). Search says 10 significant records, and the 10 in the table are consistent. However, the closing line 'shows 20 of the 38 records I pulled' contradicts the table showing only 10, and the 'baseline heat-tolerance marker' entry is described without noting it is a LI2018 baseline contrast not significant (that is fine). The 'both studies flagged for low and tank-level replication' is supported. The main issue is the inconsistent '20 shown' claim.
- ✅ **no_overclaim**: It states nothing is validated and does not imply causation, and it separates stress response from resilience.
- ✅ **criterion**: It summarizes studies, phenotypes and directions and invents no gene name or function.

**Answer** (3 queries, 17 s, $0.092)

````markdown
AREE has **38 real (non-simulated) records for gene 117690579** across 17 studies: 37 gene-expression records from 16 studies and 1 methylation record. **10 records are significant.** The clearest signal is strong up-regulation under heat stress. Overall the directions conflict: 30 up, 7 down, 1 hypo-methylated.

**Significant records (adjusted p ≤ 0.05; all 10 shown)**

| Study | Comparison | log2FC | adj. p |
|---|---|---|---|
| LI2018_HEAT_TIDAL | 35 °C 24 h vs 0 h (gill) | +6.46 | 2.2e-57 |
| GEORGE2023_HEATWAVE_PLOIDY | diploid, 30 °C + 44 °C emersion | +4.31 | 1.3e-22 |
| UW2016_HEATSHOCK_MANTLE | 40 °C 1 h (mantle, n=6) | +5.36 | 7.2e-14 |
| LI2018_HEAT_TIDAL | 35 °C 6 h | +3.89 | 1.3e-10 |
| GEORGE2023_HEATWAVE_PLOIDY | triploid, 30 °C + 44 °C emersion | +2.78 | 1.6e-08 |
| GEORGE2023_HEATWAVE_PLOIDY | diploid, 30 °C | +2.07 | 3.6e-07 |
| LIU2019_HEAT_HSF1 | 35 °C 2 h (juvenile gill) | +3.66 | 0.001 |
| LI2022_LOWSAL_VIBRIO | 10 vs 30 ppt + Vibrio, 12 h | +4.33 | 0.001 |
| DUPERRET2025_HEAT_RNASEQ | 30 vs 23 °C, spat | +1.96 | 0.038 |
| HESSER2024_VCOR | Vibrio coralliilyticus vs larvae only | −2.22 | 0.048 |

**Candidate entries (7 in total; key ones below)**
- **Heat tolerance, exposure-only, emerging** (4 studies): every study goes up, with a pooled log2FC of 4.38. But heterogeneity is very high (I² = 96%), so the size of the effect is not consistent between studies.
- **Disease resistance, high-priority cross-study, resilience-associated:** this rests on just 2 related studies, DUPERRET2025_HEAT_RNASEQ and DUPERRET2025_STARVE_RNASEQ. The starvation record alone is not significant (p = 0.087), and both studies are flagged for low and tank-level replication.
- **Disease susceptibility, emerging** (6 studies): directions are mixed (67% agree) and the pooled result is not significant (p = 0.93).
- **Baseline heat-tolerance marker** (LI2018_HEAT_TIDAL) and **methylation** (CIBNOR2021_HEAT_RRSS): one study each, neither significant.

**What this means:** the gene responds strongly and repeatedly to heat. That shows a stress response, not that animals expressing it are more resilient. Nothing here is validated, and none of it shows cause and effect.

The table above shows 20 of the 38 records I pulled. All 38 are on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## predictive-candidates, run 1: FAIL

> Are there candidates that could predict which oysters will be resilient before they are stressed?

- ✅ **facts**: ok: states 28221
- ✅ **grounded**: 1 study ids cited, all real
- ✅ **answers**: It directly answers that baseline-design candidates exist, giving the count of 28,221 and noting they are all emerging.
- ❌ **faithful**: Several claims go beyond the tool results: the effect size is labelled 'log2' though the tool only gives pooled_effect, and 'All are for thermal tolerance' is only checked for the top 15 and by zero results for other phenotypes. The study record also does not say the contrast is 'log2'.
- ✅ **no_overclaim**: It says these are associations, not causes, and that none is a validated predictor.
- ✅ **criterion**: It explains that baseline candidates exist, gives the count of 28,221, and says all are emerging and unreplicated.

**Answer** (8 queries, 22 s, $0.093)

````markdown
AREE does hold some candidates of this kind, but they are early leads, not proven predictors. 28,221 candidates come from the "baseline" design: tolerant and sensitive oysters compared *before* any stress, which is the design you'd need for a predictive marker. All of them are **emerging-tier**, meaning none has been replicated yet. All are for **thermal tolerance**. I found no baseline candidates for acidification tolerance, disease resistance or susceptibility, toxicant tolerance or recovery after stress.

**Where the top candidates come from:** All 15 top-ranked candidates come from one study, **LI2018_HEAT_TIDAL** (Li et al. 2018, *Front. Physiol.*, PRJNA407831). Its gill RNA-seq compared unstressed offspring of intertidal oysters (more heat-tolerant) with offspring of subtidal oysters (more sensitive). Here are 5 of the 15:

| Gene ID (*M. gigas*) | Effect size (log2) | Adjusted p |
|---|---|---|
| 136271467 | −9.46 | 5.7e-07 |
| 105341195 | +4.14 | 1.8e-10 |
| 105333787 | −3.14 | 6.7e-13 |
| 136271695 | +5.57 | 4.1e-07 |
| 117686755 | −1.64 | 2.7e-17 |

**Limits to keep in mind** (from the study record):
- **One study:** each candidate rests on a single study, so there's no cross-study replication. The "direction consistency = 1.0" and I² = 0 scores don't add support here, because there's only one result to compare.
- **Few groups:** only 2 populations per class (one intertidal, one subtidal at each of two sites). Those populations differ in more than tidal origin, for example growth and energy reserves.
- **Different ages:** survival was measured in 4-month-olds, but transcriptomes in 10-month-olds from the same cohorts.
- **Pooled samples:** each library pools five oysters.

These genes show a statistical association with tolerance, not a cause, and none is a validated predictor. The next step would be to test them in independent tolerant vs. sensitive families measured before stress.

The full list of 28,221 is on the [Explore page](https://sr320.github.io/AREE/explore.html) (filter contrast type = baseline_phenotype_difference), where you can also download it as a CSV.
````

## predictive-candidates, run 2: PASS

> Are there candidates that could predict which oysters will be resilient before they are stressed?

- ✅ **facts**: ok: states 28221
- ✅ **grounded**: 2 study ids cited, all real
- ✅ **answers**: It directly addresses whether baseline candidates exist for predicting resilience before stress.
- ✅ **faithful**: Claims match the tool results: 28,221 emerging candidates, LI2018_HEAT_TIDAL as the only source among the top 25, the study limitations, and WANG2020 having a baseline contrast with analysis_status not_started. Calling the effect 'log2FC' is a minor unsupported label, and the 'all for thermal tolerance' claim is inferred from the thermal_tolerance-filtered count equaling the total of 28,221.
- ✅ **no_overclaim**: It says the candidates are associations, not validated or mechanistic, and that they require replication.
- ✅ **criterion**: It says baseline candidates exist, gives the count of 28,221, and states they are all emerging and unreplicated.

**Answer** (7 queries, 21 s, $0.124)

````markdown
**Not yet any you could rely on.** AREE holds 28,221 candidates from a baseline design, where resilient and susceptible oysters were compared before any stress. That is the design a predictive marker needs. But all of them are **emerging tier**, all are for **thermal tolerance**, and none has been replicated in a second study. No baseline candidates exist for disease resistance, and none reach the high-priority cross-study tier.

**Where they come from: LI2018_HEAT_TIDAL** (Li et al. 2018, *M. gigas*, gill RNA-seq)
- The study compared offspring of intertidal oysters (more heat-tolerant) with offspring of subtidal oysters (more heat-sensitive) at 0 h, before heating. Tolerance was measured separately, as survival after 1 h at 42 °C.
- Every one of the 25 top-ranked candidates comes from this study (k_studies = 1).

**Top 5 by score (5 of 28,221):**

| Gene ID | Effect (log2FC) | Adj. p |
|---|---|---|
| 136271467 | −9.46 | 5.7e-7 |
| 105341195 | +4.14 | 1.8e-10 |
| 105333787 | −3.14 | 6.7e-13 |
| 136271695 | +5.57 | 4.1e-7 |
| 105348606 | −8.38 | 4.2e-5 |

**Limitations to keep in mind (from the study record):**
- Only two populations per group (one intertidal, one subtidal at each of two sites). The populations differ in more than tidal origin, such as growth and energy reserves.
- Each library pools five oysters.
- Survival was measured in 4-month-olds, but transcriptomes came from 10-month-olds.
- These are associations, not mechanisms or validated markers.

**What could change this:** WANG2020_TIDAL_HEAT_WGBS (methylation) also has a baseline contrast, but AREE hasn't analysed it yet. Once it is harmonized, it may be able to replicate or contradict these genes.

The full ranked list and CSV downloads are on the [Explore page](https://sr320.github.io/AREE/explore.html). Filter by contrast type "baseline_phenotype_difference".
````

## qc-failed, run 1: PASS

> Which studies failed AREE's quality control, and does their evidence count toward the results?

- ✅ **facts**: ok: mentions all of CRANDALL2022_HEAT_SEED, DUPERRET2025_POMS_STARVE, IOCAS2022_OA_ENERGY, IOCAS2024_OA_TRANSGEN, WANAMAKER2020_HEAT_SEED, WANG2024_HEAT_GILL
- ✅ **grounded**: 6 study ids cited, all real
- ✅ **answers**: It lists the failed-QC studies and says their evidence does not count toward results.
- ✅ **faithful**: The six failed studies, their reasons, n_evidence_records=0, 39% normalization disagreement and 43-48% DE all match the tool results; the claims about candidate rankings and replication counts are reasonable inferences from zero records.
- ✅ **no_overclaim**: It makes no validation or causal claims.
- ✅ **criterion**: It lists all six failed studies and states they have zero evidence records and are not in the evidence tables.

**Answer** (9 queries, 21 s, $0.234)

````markdown
**Six of AREE's 35 real studies failed quality control. None of their evidence counts.** Each has `n_evidence_records = 0`, so nothing from them appears in the evidence tables, the replication counts or the candidate rankings. AREE keeps them registered so the reason is on record and nobody curates the same dataset again.

| Study | Assay / stressor | Why it failed QC |
|---|---|---|
| CRANDALL2022_HEAT_SEED | Proteomics, temperature | One silo per arm, so the temperature effect can't be separated from the silo. The 23 °C silo also had a separate history and ciliates. |
| WANAMAKER2020_HEAT_SEED | Proteomics, temperature | No biological replication: one silo per temperature. The "replicates" are repeat runs of the same sample on the instrument. |
| DUPERRET2025_POMS_STARVE | Metabolomics, food deprivation | Fed and starved samples were run on the instrument months apart, so diet is mixed up with run batch. Two standard ways of normalizing the data disagreed on the direction for 39% of metabolites. |
| IOCAS2022_OA_ENERGY | RNA-seq, ocean acidification | A check of how much the "replicates" vary showed they aren't independent animals. 43–48% of genes came out as differentially expressed, a sign the replicates were not real. |
| IOCAS2024_OA_TRANSGEN | RNA-seq, ocean acidification | Each arm is one larval pool sequenced three times, and the two arms come from different parents, so parental pH is mixed up with genetics. |
| WANG2024_HEAT_GILL | RNA-seq, temperature | Genotype checks showed one RNA source per arm, most likely a pool, despite six separate sample labels. |

**Two other things to know:**
- The two RNA-seq studies with the same problem (IOCAS2022, WANG2024) passed AREE's first replication check, which looks for separate sample labels. The problem only showed up after the data were analyzed. Sample labels alone don't prove the samples are independent animals.
- Many studies that passed are rated "passed_with_warnings" and carry quality flags such as low replication or a single timepoint. Those flags apply to the evidence they contribute.

Several other studies are still "in_progress" or "not_started" for QC. The full study list is on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## qc-failed, run 2: FAIL

> Which studies failed AREE's quality control, and does their evidence count toward the results?

- ✅ **facts**: ok: mentions all of CRANDALL2022_HEAT_SEED, DUPERRET2025_POMS_STARVE, IOCAS2022_OA_ENERGY, IOCAS2024_OA_TRANSGEN, WANAMAKER2020_HEAT_SEED, WANG2024_HEAT_GILL
- ✅ **grounded**: 10 study ids cited, all real
- ✅ **answers**: It lists the six studies with qc_status failed and states their evidence does not count toward results.
- ❌ **faithful**: The unfinished-studies bullet says 'Every one I checked showed n_evidence_records: 0 above' for studies like VALDIVIESO2025_POMS_RS_EMSEQ and WANG2020_TIDAL_HEAT_WGBS, but only the six failed studies were checked with get_study, so this claim is unsupported.
- ✅ **no_overclaim**: No claim of validation or causation is made.
- ✅ **criterion**: It lists all six failed studies and says they contribute no evidence records (n_evidence_records: 0).

**Answer** (9 queries, 29 s, $0.238)

````markdown
**Six of AREE's 35 real studies failed quality control. None of their evidence counts.** Each one has `n_evidence_records: 0`, so nothing from them reaches the evidence tables, the cross-study replication counts or the candidate scores. AREE keeps them registered so the reason for failing is on record and the dataset isn't curated again.

| Study | Assay / question | Why it failed QC |
|---|---|---|
| CRANDALL2022_HEAT_SEED | Proteomics, spat at 29 vs 23 °C | One silo (tank) per temperature, so temperature can't be separated from silo. The 23 °C silo also had a different history and a likely ciliate outbreak |
| WANAMAKER2020_HEAT_SEED | Proteomics, spat at 29 vs 23 °C | One silo per temperature, no biological replication. Development time also dominates the proteome |
| DUPERRET2025_POMS_STARVE | Metabolomics, starved vs fed | Diet is fully confounded with the mass-spec run. Two normalization methods disagreed on the direction for 39% of metabolites |
| IOCAS2022_OA_ENERGY | RNA-seq, ocean acidification | The "replicates" weren't independent animals (effective n = 1 per group). There was no reproducible acidification signal |
| IOCAS2024_OA_TRANSGEN | RNA-seq, offspring of acidification-exposed parents | Each arm was one larval pool sequenced three times. Parental pH is confounded with which parents were used |
| WANG2024_HEAT_GILL | RNA-seq, gill after 12 h heat | A genotype check showed one RNA source per arm, so the replicates are only repeat libraries |

Things to keep in mind:
- **Studies labelled as replicated can still fail.** In IOCAS2022_OA_ENERGY and WANG2024_HEAT_GILL, the replicates looked distinct in the public metadata. The problem only showed up in AREE's checks after the reads were processed.
- **"Passed with warnings" is not the same as clean.** Most harmonized studies, such as LI2018_HEAT_TIDAL and GEORGE2023_HEATWAVE_PLOIDY, have this status. Their records do count, but they carry flags like low replication or a single timepoint, which you should weigh.
- **Unfinished studies aren't failures.** Several studies have QC "not started" or "in progress", for example VALDIVIESO2025_POMS_RS_EMSEQ and WANG2020_TIDAL_HEAT_WGBS. Every one I checked showed `n_evidence_records: 0` above, but I didn't confirm record counts for these unfinished studies.

You can see each study's full QC notes on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## thermal-replicated, run 1: PASS

> Which genes are significant for thermal tolerance in more than one real study, and do any of them disagree in direction?

- ✅ **facts**: ok: states 2919; ok: states 1389
- ✅ **grounded**: 4 study ids cited, all real
- ✅ **answers**: It gives 2,919 replicated features and 1,389 conflicting, and answers which genes replicate and which disagree.
- ✅ **faithful**: Counts, gene IDs, orthogroups, directions, studies and p-values match the tool result; the 1,389/2,919 ≈ 48% is correct arithmetic, and the 'Most of these designs are heat exposures' claim is a hedged interpretation (the tool shows temperature/multi_stressor stressors).
- ✅ **no_overclaim**: It explicitly says these are statistical associations, not validated biomarkers, and shows no causation.
- ✅ **criterion**: It gives both counts (2,919 and 1,389) and names examples of consistent genes (e.g., 117690579) and conflicting genes (e.g., 105321326).

**Answer** (1 queries, 15 s, $0.080)

````markdown
**Yes, many genes replicate, and many disagree.** In real data, **2,919 features** (genes, proteins, methylation regions or metabolites) reach adjusted p ≤ 0.05 for `thermal_tolerance` in at least 2 studies. **1,389 of them (about 48%) change in opposite directions** in different studies. I'm showing a subset of the top 25 rows below. For the full list, see the [Explore page](https://sr320.github.io/AREE/explore.html).

**Consistent across all 4 studies** (GEORGE2023_HEATWAVE_PLOIDY, LI2018_HEAT_TIDAL, LIU2019_HEAT_HSF1, UW2016_HEATSHOCK_MANTLE):

| Gene | Orthogroup | Direction | Smallest adj. p |
|---|---|---|---|
| 117690579 | none | up (7 records) | 2.2e-57 |
| 105348305 | OG0000560 | up (6) | 8.6e-56 |
| 105348304 | OG0000560 | up (6) | 1.6e-52 |
| 117680673 | OG0000560 | up (6) | 5.5e-49 |
| 105334580 | OG0020194 | up (6) | 6.4e-45 |
| 117692944 | OG0001676 | down (4) | 0.0056 |

Several other genes in orthogroup OG0000560 (105334510, 105322484, 105334234) also go up in all 4 studies.

**Conflicting directions** (all 4 studies, gill and mantle):

| Gene | Studies up / down | Smallest adj. p |
|---|---|---|
| 105321326 | 1 up / 3 down | 1.1e-54 |
| 105318591 | 1 up / 3 down | 4.2e-34 |
| 105331566 | 3 up / 1 down | 2.8e-22 |

**Things to keep in mind:**
- None of these records are simulated.
- The disagreements come from studies that differ in tissue (gill or mantle) and stressor (temperature alone, or combined with other stressors). That may explain some of the conflict, but I haven't checked study by study.
- Most of these designs are heat exposures. So these genes mostly show a **heat-stress response**, not necessarily a marker of a more tolerant animal.
- Studies are counted by study_id. Two studies from the same experiment aren't fully independent replicates.
- These are statistical associations. None is a validated biomarker, and none shows cause and effect.

I can pull up the evidence card for any of these genes, or show which ones have baseline resilient-vs-susceptible contrasts.
````

## thermal-replicated, run 2: FAIL

> Which genes are significant for thermal tolerance in more than one real study, and do any of them disagree in direction?

- ✅ **facts**: ok: states 2919; ok: states 1389
- ✅ **grounded**: 4 study ids cited, all real
- ✅ **answers**: It gives the count of replicated features, the count that conflict, and example genes of each kind.
- ❌ **faithful**: Several claims are unsupported. 'Related gene copies, not separate signals' for the shared orthogroup is an inference the tools did not give. The claim that 1,389 genes 'go up in at least one study and down in another' is not directly supported, since the conflict flag is given without per-gene details for all of them. The 'Consistent across all four studies' heading is also questionable for 117692944 (down, 4 studies), though that one does check out.
- ✅ **no_overclaim**: It states explicitly that these are associations, not causal, and that no gene is a validated biomarker.
- ✅ **criterion**: It gives 2,919 and 1,389 and names consistent examples (e.g., 117690579) and conflicting examples (e.g., 105321326).

**Answer** (1 queries, 13 s, $0.060)

````markdown
**2,919 genes** have adjusted p ≤ 0.05 for thermal tolerance in at least two real studies. **1,389 of them (about 48%) disagree in direction**: they go up in at least one study and down in another. All records are real data (none simulated). I'm showing only the top 25 rows below.

**Consistent across all four studies** (GEORGE2023_HEATWAVE_PLOIDY, LI2018_HEAT_TIDAL, LIU2019_HEAT_HSF1, UW2016_HEATSHOCK_MANTLE; gill and mantle):

| GeneID | Orthogroup | Direction (studies) | Smallest adj. p |
|---|---|---|---|
| 117690579 | – | up (4) | 2.2e-57 |
| 105348305 | OG0000560 | up (4) | 8.6e-56 |
| 105348304 | OG0000560 | up (4) | 1.6e-52 |
| 117680673 | OG0000560 | up (4) | 5.5e-49 |
| 105334580 | OG0020194 | up (4) | 6.4e-45 |
| 117692944 | OG0001676 | down (4) | 0.0056 |

Several of the strongest "up" genes share orthogroup OG0000560 (also 105334510, 105322484, 105334234). That means they're related gene copies, not separate signals.

**Examples that disagree in direction** (same four studies):

| GeneID | Studies up / down | Records up / down | Smallest adj. p |
|---|---|---|---|
| 105321326 | 1 / 3 | 2 / 4 | 1.1e-54 |
| 105318591 | 1 / 3 | 1 / 4 | 4.2e-34 |
| 105331566 | 3 / 1 | 4 / 1 | 2.8e-22 |

**Things to keep in mind:**
- Direction conflicts could come from differences in tissue (gill vs mantle), stressor (temperature vs multi-stressor) or study design. Averaging them away would hide that, so they're listed separately here.
- Studies are counted by study_id. Two studies from one experiment (for example, sharing a control group) aren't fully independent. `get_study` shows each study's design.
- These are statistical links between genes and thermal tolerance in individual studies. They don't show cause, and none of these genes is a validated biomarker.

For the full list and CSV downloads, see the [Explore page](https://sr320.github.io/AREE/explore.html). I can also pull the per-study records for any of these genes.
````

## acidification-studies, run 1: PASS

> Which studies measured acidification tolerance?

- ✅ **facts**: ok: mentions all of HARNEY2016_OA_WARM_LARVAE, IOCAS2022_OA_ENERGY, IOCAS2024_OA_TRANSGEN, LUTIER2022_OA_TIPPING, TIMMINS2014_OA_GILL
- ✅ **grounded**: 5 study ids cited, all real
- ✅ **answers**: It lists the 5 studies with acidification_tolerance, which directly answers the question.
- ✅ **faithful**: Study ids, assays, accessions, tissues, life stages, QC and analysis statuses, and the exposure_response contrast all match the tool results; 'across generations' follows from the study title, and the harmonized-study inference is hedged.
- ✅ **no_overclaim**: It says the studies show responses to exposure rather than tolerance markers and makes no validation or causal claims.
- ✅ **criterion**: The criterion is none, so it passes.

**Answer** (4 queries, 16 s, $0.061)

````markdown
AREE has **5 real studies** registered under the acidification tolerance phenotype (`acidification_tolerance`). All are on Pacific oyster (*Magallana gigas*). Searching by the ocean acidification stressor (`ocean_acidification`) finds the same 5 studies.

| Study | Assay | Tissue / life stage | Stressor(s) | QC / analysis status |
|---|---|---|---|---|
| `LUTIER2022_OA_TIPPING` | RNA-seq (PRJNA735889) | Whole animal / juvenile | Ocean acidification | Passed QC with warnings; harmonized |
| `TIMMINS2014_OA_GILL` | Proteomics (PRIDE PXD000835) | Gill / adult | Acidification + handling stress | Passed QC with warnings; harmonized |
| `HARNEY2016_OA_WARM_LARVAE` | Proteomics (PRIDE PXD002316) | Whole animal / larva | Acidification, temperature and both combined | QC and analysis not started |
| `IOCAS2022_OA_ENERGY` | RNA-seq (PRJNA826964) | Digestive gland / adult | Ocean acidification | **Failed QC**; analysis in progress |
| `IOCAS2024_OA_TRANSGEN` | RNA-seq (PRJNA1196326) | Whole larva / larva | Ocean acidification (across generations) | **Failed QC**; analysis in progress |

**How to read these:**
- **They show responses to exposure, not tolerance.** Every study compares exposed and control animals. None compares tolerant with susceptible animals, so their results show how oysters respond to acidification. They do not identify markers of a more tolerant oyster.
- **Most are labelled exposure-only.** `HARNEY2016_OA_WARM_LARVAE` is labelled resilience, and `TIMMINS2014_OA_GILL` has both labels.
- **Only two are harmonized.** `LUTIER2022_OA_TIPPING` and `TIMMINS2014_OA_GILL` are the ones most likely to have usable evidence records. I haven't checked record counts for any of the five.
- **The two IOCAS datasets** have no linked publication yet, and both failed QC.

You can see full study records and download CSVs on the [Explore page](https://sr320.github.io/AREE/explore.html). I can also pull the significant genes or proteins from the two harmonized studies.
````

## acidification-studies, run 2: PASS

> Which studies measured acidification tolerance?

- ✅ **facts**: ok: mentions all of HARNEY2016_OA_WARM_LARVAE, IOCAS2022_OA_ENERGY, IOCAS2024_OA_TRANSGEN, LUTIER2022_OA_TIPPING, TIMMINS2014_OA_GILL
- ✅ **grounded**: 5 study ids cited, all real
- ✅ **answers**: It lists the five studies that measured acidification_tolerance and describes their status and evidence.
- ✅ **faithful**: Claims match the tool results: 32,064 minus 639 gives about 31,425 for Lutier, which follows because n_studies is 2 and only those two studies have records; the directions, statuses, assays and the gene example are all supported. 'Most of the Lutier records are exposure only' is supported because that study is classified exposure_only.
- ✅ **no_overclaim**: It states that these are statistical associations and not validated biomarkers, and it does not imply causation.
- ✅ **criterion**: The criterion is 'none', so it passes.

**Answer** (7 queries, 22 s, $0.084)

````markdown
AREE has **five real studies** registered for the acidification tolerance phenotype (`acidification_tolerance`), all in Pacific oyster (*Magallana gigas*). Only **two of them have harmonized evidence so far**: 32,064 records in total.

| Study | Assay | Tissue / life stage | Status | Evidence records |
|---|---|---|---|---|
| `LUTIER2022_OA_TIPPING` (Lutier et al. 2022, *Global Change Biology*) | RNA-seq | Whole animal, juvenile | Harmonized; quality checks passed with warnings | about 31,425 |
| `TIMMINS2014_OA_GILL` (Timmins-Schiffman et al. 2014, *BMC Genomics*) | Proteomics | Gill, adult | Harmonized; quality checks passed with warnings | 639 |
| `HARNEY2016_OA_WARM_LARVAE` (Harney et al. 2016, *J Proteomics*) | Proteomics | Whole animal, larva | Not started | none yet |
| `IOCAS2022_OA_ENERGY` (BioProject PRJNA826964, no publication found) | RNA-seq | Digestive gland, adult | Quality checks failed, analysis in progress | none yet |
| `IOCAS2024_OA_TRANSGEN` (BioProject PRJNA1196326, no publication found) | RNA-seq | Whole larva | Quality checks failed, analysis in progress | none yet |

**What to keep in mind:**
- **All five compare exposed with unexposed animals.** None compares tolerant and sensitive animals, so the records show how oysters respond to low pH or high CO2. They don't show which animals are more tolerant.
- **Most of the Lutier records are labelled "exposure only".** Timmins-Schiffman also covers recovery after handling stress (`recovery_following_stress`), and some of its records are labelled "resilience".
- **Effects go in both directions**, roughly evenly up and down (16,233 up vs 15,831 down). For example, gene 105344630 is up about 1.8 log2-fold at pH ≤ 6.7 (`LUTIER2022_OA_TIPPING__ph_le6_7_vs_ph_ge7_6_23d__105344630__DESeq2`).
- **The Harney study also adds warming** on top of acidification.

These are statistical associations from individual studies, not validated biomarkers. Full tables and CSV downloads are on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## atlantic-salmon, run 1: PASS

> What does AREE show about thermal tolerance in Atlantic salmon?

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 10 study ids cited, all real
- ✅ **answers**: It says plainly that AREE has no Atlantic salmon evidence and names the species it does cover, which is what was asked.
- ✅ **faithful**: The claims match the tool results: no salmon studies, the three covered species, 10 thermal tolerance studies, and the table's assays, tissues, classes and statuses. 'Gill (juvenile)' for LIU2019 matches the registered life stage. The 'including baseline contrast' note for LI2018 also matches.
- ✅ **no_overclaim**: It states that the records are statistical associations and not validated biomarkers, and it makes no causal claims.
- ✅ **criterion**: It says AREE has no salmon evidence, names Pacific oyster and Manila clam, and keeps the oyster and clam results labelled as such rather than presenting them as salmon results.

**Answer** (3 queries, 18 s, $0.050)

````markdown
AREE has no evidence on Atlantic salmon at all, for thermal tolerance or any other phenotype. A study search for "salmon" found nothing. The only species with real evidence in this build are Pacific oyster (*Crassostrea* / *Magallana gigas*) and Manila clam (*Ruditapes philippinarum*).

AREE does hold real thermal tolerance studies in those two bivalves. There are 10 registered:

| Study | Species | Assay | Tissue | Evidence class |
|---|---|---|---|---|
| CIBNOR2021_HEAT_RRSS | *M. gigas* | Methylation | Gill | Resilience |
| GEORGE2023_HEATWAVE_PLOIDY | *M. gigas* | RNA-seq | Gill | Stress response |
| LI2018_HEAT_TIDAL | *M. gigas* | RNA-seq | Gill | Resilience + stress response (includes a baseline resilient-vs-susceptible contrast) |
| LIU2019_HEAT_HSF1 | *M. gigas* | RNA-seq | Gill (juvenile) | Stress response |
| UW2016_HEATSHOCK_MANTLE | *M. gigas* | RNA-seq | Mantle | Exposure only |
| WANG2020_TIDAL_HEAT_WGBS | *M. gigas* | Methylation | Gill | Resilience + stress response (analysis not started) |
| WANG2024_HEAT_GILL | *M. gigas* | RNA-seq | Gill | Exposure only (QC failed, in progress) |
| WANG2026_ERK_HEAT | *M. gigas* | Proteomics | Gill | Exposure only + resilience (analysis not started) |
| PERUZZA2023_HEATWAVE | *R. philippinarum* | RNA-seq | Digestive gland | Exposure only |
| TUCCI2025_HEAT_PRIMING | *R. philippinarum* | RNA-seq | Digestive gland | Exposure only + resilience |

Even in oysters, these are statistical associations from individual studies, not validated biomarkers.

If oyster or clam results would still help, I can look up:
- Genes that are significant for thermal tolerance in more than one study, and whether their directions conflict
- Top-ranked thermal tolerance candidates, such as resilience-associated or baseline-contrast candidates

You can see the full study tables on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## atlantic-salmon, run 2: PASS

> What does AREE show about thermal tolerance in Atlantic salmon?

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 10 study ids cited, all real
- ✅ **answers**: It says plainly that AREE holds no Atlantic salmon evidence and lists the species it does cover.
- ✅ **faithful**: The claims match the tool results: 35 real studies, 0 salmon studies, three species, 10 thermal tolerance studies, and the per-study details (QC failed for WANG2024, not_started for WANG2020 and WANG2026). The descriptions such as 'Failed quality control; analysis in progress' and 'Not analysed yet' follow from the qc and analysis statuses.
- ✅ **no_overclaim**: It makes no validation or causal claims.
- ✅ **criterion**: It states AREE has no salmon evidence, names Pacific oyster and Manila clam, and labels the oyster and clam studies as comparison material rather than salmon results.

**Answer** (3 queries, 16 s, $0.049)

````markdown
AREE has no evidence on Atlantic salmon. None of its 35 real studies are on salmon, and a species search for "salmon" finds no studies. All of its current evidence comes from three bivalve species: Pacific oyster (*Crassostrea gigas* / *Magallana gigas*) and Manila clam (*Ruditapes philippinarum*).

AREE does have 10 thermal tolerance studies for those bivalves, if they're useful as a comparison:

| Study | Species | Data type | What the study compares |
|---|---|---|---|
| LI2018_HEAT_TIDAL | Pacific oyster | RNA-seq, gill | Resilient vs. susceptible animals before stress, plus the response to heat |
| CIBNOR2021_HEAT_RRSS | Pacific oyster | Methylation, gill | Heat-resistant vs. heat-susceptible families during heat exposure |
| GEORGE2023_HEATWAVE_PLOIDY | Pacific oyster | RNA-seq, gill | Response to heatwaves and combined stressors |
| LIU2019_HEAT_HSF1 | Pacific oyster | RNA-seq, gill | Response to heat |
| UW2016_HEATSHOCK_MANTLE | Pacific oyster | RNA-seq, mantle | Exposure to heat only |
| WANG2020_TIDAL_HEAT_WGBS | Pacific oyster | Methylation, gill | Not analysed yet |
| WANG2024_HEAT_GILL | Pacific oyster | RNA-seq, gill | Failed quality control; analysis in progress |
| WANG2026_ERK_HEAT | Pacific oyster | Proteomics, gill | Not analysed yet |
| PERUZZA2023_HEATWAVE | Manila clam | RNA-seq, digestive gland | Exposure to heat only |
| TUCCI2025_HEAT_PRIMING | Manila clam | RNA-seq, digestive gland | Heat exposure, and how primed vs. unprimed clams differ under it |

Other phenotypes with evidence in AREE are acidification tolerance, disease resistance and susceptibility, recovery after stress, and toxicant tolerance.

I can't tell you which genes are linked to thermal tolerance in salmon, because AREE has no salmon data. For salmon studies you'd need to look at other sources.

Some questions I can answer:
- Which genes are significant for thermal tolerance in more than one oyster study, and do any of them change in opposite directions?
- What are the top thermal tolerance candidate markers from studies comparing resilient and susceptible oysters before stress?

You can also browse full tables on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## follow-up-raw-data, run 1: PASS

> Which studies measured acidification tolerance?
> Of those, which have raw data available for reanalysis?

- ✅ **facts**: ok: mentions all of IOCAS2022_OA_ENERGY, IOCAS2024_OA_TRANSGEN, LUTIER2022_OA_TIPPING, TIMMINS2014_OA_GILL
- ✅ **grounded**: 5 study ids cited, all real
- ✅ **answers**: It identifies which of the five acidification studies have raw data (four) and which do not (HARNEY2016), as asked.
- ✅ **faithful**: Raw/processed flags, accessions, analysis modes and QC failed/in progress all match the tool results; 'raw reads are public' is supported by the citation text (reads public dates), and the 'only say failed' remark is consistent with the records.
- ✅ **no_overclaim**: No claims of validation or causation; the answer is hedged about what reanalysis could show.
- ✅ **criterion**: It covers the same five acidification studies from the prior answer and states which have raw data and which do not.

**Answer** (4 queries, 23 s, $0.100)

````markdown
**Four of the five studies have raw data available.** The exception is `HARNEY2016_OA_WARM_LARVAE`, which has processed results only.

| Study | Raw data | Processed data | Where the raw data is | Analysis mode |
|---|---|---|---|---|
| `LUTIER2022_OA_TIPPING` | ✅ | ❌ | NCBI BioProject PRJNA735889 (RNA-seq) | raw reanalysis |
| `TIMMINS2014_OA_GILL` | ✅ | ✅ | PRIDE PXD000835 (proteomics) | raw reanalysis |
| `IOCAS2022_OA_ENERGY` | ✅ | ❌ | NCBI BioProject PRJNA826964 (RNA-seq) | raw reanalysis |
| `IOCAS2024_OA_TRANSGEN` | ✅ | ❌ | NCBI BioProject PRJNA1196326 (RNA-seq) | raw reanalysis |
| `HARNEY2016_OA_WARM_LARVAE` | ❌ | ✅ | PRIDE PXD002316 is listed, but AREE marks raw data unavailable | processed results only |

**If you plan to reanalyse these:**
- **The two IOCAS datasets failed QC in AREE's pipeline** and are still in progress. The raw reads are public, but they may need extra checking before use. Their study records only say "failed", not why.
- **`TIMMINS2014_OA_GILL` is the only study with both raw and processed data.** That lets you check a reanalysis against the authors' own results.
- **All four still compare exposed with unexposed animals.** Reanalysing them can show responses to acidification, but not differences between tolerant and susceptible animals, unless the original designs hold information on individual phenotypes. AREE's records don't show that they do.

I can open the full registration record for any of these, with sample sizes, comparisons and listed limitations, to help you judge whether it's worth reanalysing.
````

## follow-up-raw-data, run 2: PASS

> Which studies measured acidification tolerance?
> Of those, which have raw data available for reanalysis?

- ✅ **facts**: ok: mentions all of IOCAS2022_OA_ENERGY, IOCAS2024_OA_TRANSGEN, LUTIER2022_OA_TIPPING, TIMMINS2014_OA_GILL
- ✅ **grounded**: 5 study ids cited, all real
- ✅ **answers**: It states which of the five acidification studies have raw data available and which do not, as the question and criterion require.
- ✅ **faithful**: Raw-data flags, accessions, QC and analysis statuses all match the tool results; the claim that the IOCAS reads are public follows from the citation text. The advice to start with read quality is general guidance, and AREE does not present it as a result.
- ✅ **no_overclaim**: It makes no claims of validation, causation or mechanism.
- ✅ **criterion**: It refers to the five studies from the previous answer and splits them into four with raw data and HARNEY2016_OA_WARM_LARVAE without.

**Answer** (5 queries, 23 s, $0.085)

````markdown
Four of the five studies have raw data available. Only HARNEY2016_OA_WARM_LARVAE doesn't.

| study_id | Raw data | Where the data are | AREE status |
|---|---|---|---|
| LUTIER2022_OA_TIPPING | Yes | RNA-seq, NCBI BioProject PRJNA735889 | Reanalyzed from raw data and harmonized. QC passed with warnings |
| TIMMINS2014_OA_GILL | Yes (processed results also available) | Proteomics, PRIDE PXD000835 | Reanalyzed from raw data and harmonized. QC passed with warnings |
| IOCAS2022_OA_ENERGY | Yes | RNA-seq, NCBI BioProject PRJNA826964 | Raw reanalysis in progress. **QC failed** |
| IOCAS2024_OA_TRANSGEN | Yes | RNA-seq, NCBI BioProject PRJNA1196326 | Raw reanalysis in progress. **QC failed** |
| HARNEY2016_OA_WARM_LARVAE | **No** (processed results only) | Proteomics, PRIDE PXD002316 | Harmonized from published processed results, not raw data. QC not started |

**Things to keep in mind:**
- **The two IOCAS datasets failed AREE's QC.** The raw reads are public, but they didn't pass quality checks. If you reanalyze them, start by looking at read quality. AREE's QC details are in each study's registration record.
- **Only two of the four have finished reanalysis:** LUTIER2022_OA_TIPPING and TIMMINS2014_OA_GILL.
- **AREE doesn't list raw data for HARNEY2016_OA_WARM_LARVAE.** It has a PRIDE accession, but AREE uses only its processed results.

I can pull up the full registration record for any of these, including sample sizes, limitations and QC notes for each comparison.
````

## li2018-design, run 1: PASS

> What did the LI2018_HEAT_TIDAL study compare, and in which tissue?

- ✅ **facts**: ok: mentions all of gill
- ✅ **grounded**: 1 study ids cited, all real
- ✅ **answers**: It states the three registered comparisons with treatment and control conditions and identifies gill as the tissue.
- ✅ **faithful**: All claims match the study record, including the comparisons, n values, 87,698 records, the LTJ02 exclusion, the limitations, the QC status and the flags; 'more heat tolerant' and the 'four wild populations' design are supported by the record.
- ✅ **no_overclaim**: It says the tolerance contrast shows an association, not a mechanism, and makes no validation or causal claims.
- ✅ **criterion**: It describes each comparison's treatment and control conditions and the gill tissue, consistent with the registration.

**Answer** (1 queries, 12 s, $0.055)

````markdown
**LI2018_HEAT_TIDAL** used RNA-seq to compare gene expression in the **gill** of adult (10-month-old) *Magallana gigas*. The study has three comparisons, and all of them use gill tissue. It is real data (simulated=false), based on Li et al. 2018, *Front. Physiol.* 9:825 (BioProject PRJNA407831).

| Comparison | What was compared | Evidence class |
|---|---|---|
| `heat_35c_24h_vs_0h` | 24 h in 35 °C seawater vs. the same tanks before heating (0 h), n=23. This is the primary comparison for meta-analysis. | Stress response (exposure response) |
| `heat_35c_6h_vs_0h` | 6 h in 35 °C seawater vs. 0 h, n=23 | Stress response (exposure response) |
| `intertidal_vs_subtidal_0h` | Offspring of intertidal populations (more heat tolerant) vs. offspring of subtidal populations, before any heating, n=11 | Resilience-associated (baseline phenotype difference) |

**Design:** The oysters were first-generation offspring of four wild Bohai Sea populations, one intertidal and one subtidal at each of two sites (Bayuquan and Laoting). All were raised in a common garden. Each RNA library pooled five oysters. The authors excluded one outlier sample (LTJ02), and AREE followed that exclusion. In total the study contributed 87,698 evidence records.

**Limitations to keep in mind** (from the registration):
- The heat comparisons have no unheated group sampled at the same time, so the effect of heat can't be separated from time and handling.
- The tolerance contrast rests on only two populations per tidal class, and those populations differ in more than tidal origin. It shows an association, not a mechanism.
- Heat tolerance (survival after 1 h at 42 °C) was measured in 4-month-olds, but the transcriptomes come from 10-month-olds of the same cohorts.
- QC status is "passed with warnings" (single timepoint, replication at the family level).

The full comparison tables are on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## li2018-design, run 2: PASS

> What did the LI2018_HEAT_TIDAL study compare, and in which tissue?

- ✅ **facts**: ok: mentions all of gill
- ✅ **grounded**: 1 study ids cited, all real
- ✅ **answers**: It states the three comparisons with their treatment and control conditions and names gill as the tissue.
- ✅ **faithful**: All claims match the study record: gill tissue, 10-month-old F1 oysters, the three comparisons, n values (23, 23, 11), survival after 42 °C for 1 h, the limitations, 87,698 records, and the accessions. Calling the oysters 'adult' follows the registered life_stage. The 'libraries' wording for sample_size is a reasonable reading of the record.
- ✅ **no_overclaim**: It says the tolerance comparison shows an association with resilience, not a mechanism, and makes no validation or causal claims.
- ✅ **criterion**: It describes each comparison's treatment and control conditions as registered, and the tissue (gill), consistent with the study record.

**Answer** (1 queries, 10 s, $0.038)

````markdown
**LI2018_HEAT_TIDAL** used RNA-seq on **gill tissue** from adult Pacific oysters (*Magallana gigas*), 10 months old. The oysters were the first generation (F1) bred from intertidal and subtidal wild populations at two sites in the Bohai Sea, China, and raised together in the same conditions. AREE re-analysed the raw data, and it holds three comparisons:

| Comparison | What was compared | Contrast type | Evidence class | n |
|---|---|---|---|---|
| `heat_35c_24h_vs_0h` (main comparison for meta-analysis) | 24 h in 35 °C seawater vs the same tanks before heating | exposure_response | stress_response | 23 libraries |
| `heat_35c_6h_vs_0h` | 6 h in 35 °C vs before heating | exposure_response | stress_response | 23 libraries |
| `intertidal_vs_subtidal_0h` | Before any heating: offspring of intertidal (heat-tolerant) vs subtidal populations | baseline_phenotype_difference | resilience | 11 libraries, 2 populations per group |

All three are mapped to stressor **temperature** and phenotype **thermal_tolerance**. Tolerance was measured as survival after 1 h at 42 °C, which was higher in the intertidal offspring at both sites.

**Limitations to keep in mind** (from the study record):
- **No unheated control group:** the heat comparisons use the same tanks before heating as the control, so time and handling can't be separated from the heat effect.
- **Pooled samples:** each library pools five oysters.
- **Only two populations per group in the tolerance comparison:** the populations differ in more than tidal origin, so this is evidence of an association with resilience, not of a mechanism.
- **Different ages:** survival was measured at 4 months, gene expression at 10 months.
- **Data source:** quality check passed with warnings. The authors' outlier sample (LTJ02) was left out.

The study contributes **87,698 evidence records** across the three comparisons (DOI 10.3389/fphys.2018.00825; BioProject PRJNA407831). It's a real study, not simulated data. You can browse its records on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## low-salinity, run 1: PASS

> Which studies looked at low salinity, and what phenotype did they measure?

- ✅ **facts**: ok: mentions all of LI2022_LOWSAL_VIBRIO
- ✅ **grounded**: 1 study ids cited, all real
- ✅ **answers**: It names the single low-salinity study, LI2022_LOWSAL_VIBRIO, and gives its registered phenotype, disease_susceptibility.
- ✅ **faithful**: All claims match the tool results: 10 vs 30 ppt, 12 h and 48 h, 100% vs 33% mortality, 51,197 records, the limitations, the QC status and the quality flag. The claim that no multi-stressor study includes salinity is supported by the four listed studies' stressors. The 'disease-associated' evidence class follows from the resilience_classification of 'disease'.
- ✅ **no_overclaim**: It describes the records as statistical associations and explicitly says they are not validated markers; no causal claims are made.
- ✅ **criterion**: It names LI2022_LOWSAL_VIBRIO and its registered phenotype, disease_susceptibility, for both comparisons.

**Answer** (6 queries, 16 s, $0.069)

````markdown
Only one real study in AREE looked at low salinity: **LI2022_LOWSAL_VIBRIO**. It measured **disease susceptibility**, recorded as mortality after a *Vibrio* infection. No real studies are registered under "freshwater exposure", and none of the four multi-stressor studies includes salinity.

| | LI2022_LOWSAL_VIBRIO |
|---|---|
| Species / tissue / stage | *Magallana gigas*, gill, adult |
| Assay | RNA-seq, reanalysed by AREE from the raw reads (BioProject PRJNA756710) |
| Design | All oysters were injected with *Vibrio alginolyticus*, then held at 10 ppt (low salinity) or 30 ppt (normal) |
| Comparisons | 10 vs 30 ppt at 12 h and at 48 h after infection |
| Phenotype | Disease susceptibility (it goes up at low salinity). Over 12 days, 100% of the 10 ppt group died, against 33% of the 30 ppt group |
| Evidence class | Disease-associated, exposure response |
| Evidence records | 51,197 |
| Citation | Li et al. 2022, *Frontiers in Immunology* 13:859975 (doi:10.3389/fimmu.2022.859975) |

**Things to keep in mind:**
- **Salinity was never tested on its own.** Every oyster was infected, and no uninfected oysters were sequenced. The comparisons show how low salinity changes gene expression in *infected* oysters, not the effect of salinity alone.
- **Mortality was counted per group, not per sequenced oyster.** The oysters that were sequenced were sampled at 12 h and 48 h, and nobody knows whether they would have died.
- **The sample size is small.** Each comparison is 3 tanks against 3 tanks, with one pooled library of 3 oysters per tank. The same tanks were sampled at both time points, so the 12 h and 48 h results are not independent.
- The study passed quality control with warnings and is flagged for an ambiguous phenotype definition.

These records are statistical associations from one study, not validated markers of salinity tolerance. The full evidence tables are on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## low-salinity, run 2: PASS

> Which studies looked at low salinity, and what phenotype did they measure?

- ✅ **facts**: ok: mentions all of LI2022_LOWSAL_VIBRIO
- ✅ **grounded**: 1 study ids cited, all real
- ✅ **answers**: It names the single low-salinity study, LI2022_LOWSAL_VIBRIO, and its registered phenotype, disease_susceptibility.
- ✅ **faithful**: Every claim traces to the tool results: 51,197 records, the 10 vs 30 ppt comparisons at 12 h and 48 h, 100% vs 33% mortality, the pooled 3-vs-3 tank design, the QC status and the quality flag. The gloss 'death after Vibrio infection' is a reasonable reading of the percent-mortality unit.
- ✅ **no_overclaim**: It states that nothing here is replicated and that the evidence reflects a disease or exposure response, with no claim of validation or causation.
- ✅ **criterion**: It names the study and the phenotype as registered, disease_susceptibility, for both comparisons.

**Answer** (5 queries, 18 s, $0.054)

````markdown
Low salinity may be coded as salinity or freshwater_exposure.AREE has **one** real study on low salinity: **LI2022_LOWSAL_VIBRIO**. The phenotype it measured is **disease susceptibility**, meaning death after *Vibrio* infection, not salinity tolerance. I searched both the `salinity` and `freshwater_exposure` stressor terms. No study is tagged `freshwater_exposure`.

| | LI2022_LOWSAL_VIBRIO |
|---|---|
| Species / tissue | *Magallana gigas*, adult gill |
| Assay | RNA-seq, reanalysed by AREE from the raw reads (BioProject PRJNA756710) |
| Comparisons | 10 ppt vs 30 ppt seawater, sampled 12 h and 48 h after injection with *Vibrio alginolyticus* |
| Phenotype | `disease_susceptibility` (increase) |
| Evidence class | disease, exposure_response contrast |
| Records | 51,197 harmonized evidence records |
| Citation | Li et al. 2022, *Front. Immunol.* 13:859975 |

**Caveats from the study record:**
- **Both groups were infected.** No oysters without an injection were sequenced. So the comparisons show how low salinity changes gene expression in *infected* oysters, not the effect of salinity on its own.
- **Mortality was counted per group,** not per sequenced oyster: 100% at 10 ppt vs 33% at 30 ppt over 12 days. The oysters that were sequenced were sampled at 12 h and 48 h, and what happened to them afterwards is unknown.
- **The sample is small and pooled.** Each comparison is 3 tanks vs 3 tanks, one library per tank made from 3 pooled oysters. The same tanks were sampled at both timepoints, so the two comparisons are not independent.
- **Data quality:** QC passed with warnings, and the study is flagged `ambiguous_phenotype_definition`.

Any gene from this study is evidence of a disease or exposure response, not a marker of salinity-resilient oysters. With only one study, nothing here is replicated.

If it helps, I can show the most significant genes from either timepoint. Full tables and CSV downloads are on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## metabolites, run 1: PASS

> Which metabolites are associated with thermal tolerance in oysters?

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 10 study ids cited, all real
- ✅ **answers**: It says plainly that AREE holds no metabolite evidence for thermal tolerance and names what evidence does exist.
- ✅ **faithful**: Claims match the tool results: feature types with evidence are gene/protein/methylation_region, 10 studies, assay and species groupings, QC statuses, and the LI2018 baseline contrast. The claim that no real study has produced metabolite records follows from the overview's feature_types list.
- ✅ **no_overclaim**: No gene or candidate is called validated or causal, and it notes that exposure-response results show a stress response rather than resilience.
- ✅ **criterion**: It states there are no metabolite records and names no metabolites as findings.

**Answer** (3 queries, 15 s, $0.052)

````markdown
**None at the moment.** AREE has no metabolite evidence for thermal tolerance, or for any other phenotype. Metabolite features are part of AREE's vocabulary (`metabolite_feature`), but no registered real study has produced metabolite records. The only feature types with evidence are **genes, proteins and methylation regions**.

**What AREE does hold for thermal tolerance:** 10 real studies, using the omics types below.

| Assay | Studies | Species |
|---|---|---|
| RNA-seq (genes) | GEORGE2023_HEATWAVE_PLOIDY, LI2018_HEAT_TIDAL, LIU2019_HEAT_HSF1, UW2016_HEATSHOCK_MANTLE, WANG2024_HEAT_GILL | *M. gigas* |
| RNA-seq (genes) | PERUZZA2023_HEATWAVE, TUCCI2025_HEAT_PRIMING | *R. philippinarum* (Manila clam) |
| Methylation | CIBNOR2021_HEAT_RRSS, WANG2020_TIDAL_HEAT_WGBS | *M. gigas* |
| Proteomics | WANG2026_ERK_HEAT | *M. gigas* |

Keep in mind:
- Some of these studies aren't fully processed yet. WANG2024_HEAT_GILL failed quality control and is still being analysed. WANG2020_TIDAL_HEAT_WGBS and WANG2026_ERK_HEAT haven't been analysed yet.
- The designs differ. Some studies compare resilient with susceptible animals (for example, the baseline contrast in LI2018_HEAT_TIDAL). Others only measure how animals respond to heat exposure. Exposure-response results show a stress response, not which animals are resilient.

**Questions I can answer instead:**
- Which genes or proteins are significant for thermal tolerance in more than one study, and do their directions agree?
- What are the top-ranked thermal-tolerance candidates? I can also limit this to studies comparing resilient and susceptible animals before any stress, which is the design needed for a predictive marker.

You can browse the full study tables and download CSVs on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## metabolites, run 2: PASS

> Which metabolites are associated with thermal tolerance in oysters?

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 10 study ids cited, all real
- ✅ **answers**: It states plainly that AREE holds no metabolite records for thermal tolerance, which the get_evidence result (0 records) supports.
- ✅ **faithful**: The claims about study counts (10 studies, 8 oyster and 2 clam), assay types, QC and analysis status, and the resilience classifications match the list_studies output; 'none used metabolomics' follows from the assay types listed for thermal tolerance and the feature types in the overview, though 'for any other phenotype' is only loosely supported by the overview's feature_types list, which is sufficient.
- ✅ **no_overclaim**: It calls nothing validated or causal and notes candidates are statistical associations to follow up.
- ✅ **criterion**: It says there is no metabolite evidence and names no metabolites as findings.

**Answer** (3 queries, 14 s, $0.050)

````markdown
AREE has **no metabolite evidence for thermal tolerance in oysters**, or for any other phenotype. A search for metabolite records under `thermal_tolerance` found 0 records. The current build only holds three kinds of features: **genes, methylation regions and proteins**. None of the registered studies used metabolomics.

**What AREE does have for thermal tolerance:** 10 real studies, 8 in *Magallana gigas* and 2 in the Manila clam *Ruditapes philippinarum*.

| Data type | Oyster studies |
|---|---|
| RNA-seq (genes) | LI2018_HEAT_TIDAL, GEORGE2023_HEATWAVE_PLOIDY, LIU2019_HEAT_HSF1, UW2016_HEATSHOCK_MANTLE, WANG2024_HEAT_GILL (QC failed; analysis still in progress) |
| DNA methylation | CIBNOR2021_HEAT_RRSS, WANG2020_TIDAL_HEAT_WGBS (analysis not started) |
| Proteomics | WANG2026_ERK_HEAT (analysis not started) |

The two clam studies (PERUZZA2023_HEATWAVE, TUCCI2025_HEAT_PRIMING) are RNA-seq.

These studies measure different things, so the categories matter:
- **Resilience-associated comparisons** (tolerant vs. susceptible animals): CIBNOR2021_HEAT_RRSS, the baseline comparison in LI2018_HEAT_TIDAL, and part of TUCCI2025_HEAT_PRIMING (clam).
- **Stress-response or exposure-only comparisons**: the others. These show how animals respond to heat, not which animals are tolerant.

**Next steps you could take:**
- Ask which **genes or methylation regions** are significant for thermal tolerance in more than one study, and whether their directions agree.
- Ask for the ranked **candidate biomarkers** for thermal tolerance. Any candidate is a statistical association to follow up, not a validated marker.
- Browse the full tables on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## multi-omics, run 1: PASS

> Are there multi-omics convergence candidates for any phenotype?

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 0 study ids cited, all real
- ✅ **answers**: It directly says there are no multi-omics convergence candidates, matching the 0-result tool output.
- ✅ **faithful**: Tier counts (4,962; 270,270), 1,236,316 records, 35 real studies, feature types and phenotypes all match the overview; the metabolite_linked count of 0 is a reasonable inference since the overview's candidate_tiers lists no such tier, and 'no metabolite records' follows from feature_types.
- ✅ **no_overclaim**: It states explicitly that no candidate is validated and that emerging-tier candidates need replication; the speculation about feature types is hedged and flagged as unchecked.
- ✅ **criterion**: It says plainly there are currently no multi-omics convergence candidates and invents none.

**Answer** (2 queries, 13 s, $0.022)

````markdown
**No. AREE's current build has no multi-omics convergence candidates for any phenotype.** When I filtered candidates to the `multi_omics_convergence` tier, I got 0 results.

**What the build does contain** (real studies only):

| Candidate tier | Count |
|---|---|
| `high_priority_cross_study` | 4,962 |
| `emerging` | 270,270 |
| `multi_omics_convergence` | 0 |
| `metabolite_linked_convergence` | 0 |

These candidates come from 1,236,316 evidence records across 35 real studies. The evidence covers three feature types: genes, methylation regions and proteins. There are no metabolite records, which explains why the metabolite-linked tier is empty. The tier counts alone don't explain why no candidate combines different feature types (for example, a gene and its methylation region, or a transcript and its protein). One possibility is that the feature types rarely come from the same comparisons or contexts. I haven't checked that.

**Phenotypes with evidence:** acidification_tolerance, disease_resistance, disease_susceptibility, recovery_following_stress, thermal_tolerance and toxicant_tolerance.

**Some next steps I can help with:**
- List the `high_priority_cross_study` candidates for one phenotype, such as thermal_tolerance. These are features replicated across studies.
- Check one gene for both expression and methylation records, and whether their directions agree. This would be an informal look, not a formal convergence call.
- Find features significant in several studies and flag any that conflict in direction.

None of these candidates is a validated biomarker, whatever its tier or score. Emerging-tier candidates in particular still need replication. For full candidate tables and CSV downloads, see the Explore page: https://sr320.github.io/AREE/explore.html
````

## multi-omics, run 2: PASS

> Are there multi-omics convergence candidates for any phenotype?

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 0 study ids cited, all real
- ✅ **answers**: It plainly states there are no multi-omics convergence candidates, directly answering the question.
- ✅ **faithful**: Counts (4,962, 270,270, 0), 35 real studies, feature types, phenotypes and species all match tool results; the metabolite-linked 0 follows from tier_counts omitting it, and the absence of metabolite records follows from the feature_types list.
- ✅ **no_overclaim**: It says the candidates are statistical associations, not validated biomarkers, and that emerging candidates need replication; it makes no causal claims.
- ✅ **criterion**: It states there are currently no multi-omics convergence candidates and does not invent any.

**Answer** (2 queries, 12 s, $0.017)

````markdown
No. AREE's current build has **no multi-omics convergence candidates for any phenotype**. A search of that tier came back empty, and it doesn't appear among the build's candidate tiers at all.

**Candidates in the current build (real studies only)**

| Tier | Candidates |
|---|---|
| high_priority_cross_study | 4,962 |
| emerging | 270,270 |
| multi_omics_convergence | 0 |
| metabolite_linked_convergence | 0 |

**Why there are none**

The evidence comes from 35 real studies and covers three feature types: genes, proteins and methylation regions. There are no metabolite records. That explains why the metabolite-linked tier is empty. Gene, protein and methylation data do exist, but no feature currently has the cross-layer support needed to qualify as multi-omics convergence. The tools don't report the exact qualifying rule, so I can't say which specific condition fails.

**What does have evidence**

- **Phenotypes:** acidification_tolerance, disease_resistance, disease_susceptibility, recovery_following_stress, thermal_tolerance, toxicant_tolerance
- **Species:** *Crassostrea / Magallana gigas*, *Ruditapes philippinarum*

**A caution on the other tiers**

- **High priority cross-study:** these are the strongest candidates, with replication across studies. They are still statistical associations, not validated biomarkers.
- **Emerging:** these need replication before you rely on them.

**Questions I can answer next**

- The top high priority cross-study candidates for one phenotype, such as thermal_tolerance.
- Whether a gene shows up at more than one omics layer. I can compare its gene, protein and methylation-region records, though AREE doesn't treat that as a formal convergence score.

You can browse full candidate tables and download CSVs on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## off-topic-code, run 1: PASS

> Write me a Python function that sorts a list of dictionaries by a key.

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 0 study ids cited, all real
- ✅ **answers**: The request is off-topic and the answer declines it, which is the expected handling.
- ✅ **faithful**: No tools were called and the answer makes no claims about AREE's data; the suggested question is only a prompt.
- ✅ **no_overclaim**: The answer does not call any gene or candidate validated, proven or causal.
- ✅ **criterion**: It declines briefly in one sentence, says it only answers questions about AREE's evidence, does not write the code, and suggests an answerable question.

**Answer** (0 queries, 2 s, $0.003)

````markdown
I can only answer questions about AREE's evidence on aquaculture resilience, so I can't help with Python code. Here's a question I could answer instead: *"Which genes are significant for thermal tolerance in more than one Pacific oyster study, and do any of them change in opposite directions?"*
````

## off-topic-code, run 2: PASS

> Write me a Python function that sorts a list of dictionaries by a key.

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 0 study ids cited, all real
- ✅ **answers**: The request is off-topic and the answer declines it and suggests an in-scope question, as the criterion requires.
- ✅ **faithful**: No claims about AREE data are made; the suggested question is only a suggestion.
- ✅ **no_overclaim**: No genes or candidates are called validated or causal.
- ✅ **criterion**: It declines briefly, says it only answers questions about AREE's evidence, and does not write the code.

**Answer** (0 queries, 2 s, $0.002)

````markdown
I can only help with questions about AREE's evidence, so I can't write general-purpose code. You could ask me something like *"Which genes are significant for thermal tolerance in more than one study, and do any of them disagree in direction?"*
````

## off-topic-poem, run 1: PASS

> Write a short poem about oysters.

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 0 study ids cited, all real
- ✅ **answers**: The answer declines the off-topic poem request and suggests an in-scope question, as the scope instructions require.
- ✅ **faithful**: No claims about AREE's data are made; the answer only states its scope and suggests a question.
- ✅ **no_overclaim**: No genes or candidates are called validated or causal.
- ✅ **criterion**: It declines briefly, explains it only answers questions about AREE's evidence, and does not write the poem.

**Answer** (0 queries, 3 s, $0.003)

````markdown
I can only answer questions about AREE's oyster and aquaculture evidence, so I can't write poems. You could ask me something like *"Which genes are linked to thermal tolerance in Pacific oyster in more than one study?"*
````

## off-topic-poem, run 2: PASS

> Write a short poem about oysters.

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 0 study ids cited, all real
- ✅ **answers**: It declines the off-topic poem request and suggests an in-scope question, as the scope rules require.
- ✅ **faithful**: No tools were called and the answer makes no claims about AREE's data; it only suggests a question.
- ✅ **no_overclaim**: It does not call any gene or candidate validated, proven or causal.
- ✅ **criterion**: It declines briefly, says it only answers questions about AREE's evidence, and does not write the poem.

**Answer** (0 queries, 2 s, $0.002)

````markdown
I can only answer questions about AREE's evidence, so I can't write poems. I can help with oyster questions like this one, though: *"Which genes are significant for thermal tolerance in more than one Pacific oyster study, and do any of them disagree in direction?"*
````

## phenotypes, run 1: PASS

> Which phenotypes does AREE have evidence for?

- ✅ **facts**: ok: mentions all of acidification_tolerance, disease_resistance, disease_susceptibility, recovery_following_stress, thermal_tolerance, toxicant_tolerance
- ✅ **grounded**: 0 study ids cited, all real
- ✅ **answers**: It lists the six phenotypes with evidence, which directly answers the question.
- ✅ **faithful**: Counts, phenotypes, stressors, species, feature types and tiers all match the overview result; 'Manila clam' and 'Pacific oyster' are acceptable background naming, and the per-phenotype descriptions are plain-language glosses of the ids.
- ✅ **no_overclaim**: It states that no candidate is validated, that a high score does not show cause and effect, and that emerging candidates need replication.
- ✅ **criterion**: The criterion is none, so this passes.

**Answer** (1 queries, 11 s, $0.037)

````markdown
AREE has real evidence for **six resilience-related phenotypes**. This comes from 35 real studies and about 1.24 million harmonized evidence records:

| Phenotype (controlled id) | What it covers |
|---|---|
| `thermal_tolerance` | Tolerating heat or temperature stress |
| `acidification_tolerance` | Tolerating ocean acidification |
| `disease_resistance` | Resisting infection |
| `disease_susceptibility` | Being vulnerable to infection |
| `recovery_following_stress` | Recovering after a stress event |
| `toxicant_tolerance` | Tolerating pollutants or toxins |

**Other useful context:**
- **Stressors with evidence:** temperature, ocean acidification, salinity, pathogen challenge, pollutant exposure, harmful algal bloom, nutritional limitation, handling/transport stress, and multi-stressor designs.
- **Species:** Pacific oyster (listed under both *Crassostrea gigas* and *Magallana gigas*) and Manila clam (*Ruditapes philippinarum*).
- **Feature types:** genes, proteins and DNA methylation regions.
- **Candidate biomarkers:** 275,232 in total. Of these, 4,962 are "high-priority cross-study" and 270,270 are "emerging", meaning they need replication before anyone relies on them. No candidate is a validated biomarker, and a high score doesn't show cause and effect.

These numbers leave out AREE's 7 simulated demo studies, which are excluded here.

**Possible next questions:**
- Which studies cover a given phenotype (for example, thermal tolerance)?
- Which genes are significant for a phenotype in more than one study, and do any of them disagree in direction?
- What are the top candidates for a phenotype, kept separate by evidence class (resilience-associated or stress-response)?

For full tables and CSV downloads, see the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## phenotypes, run 2: PASS

> Which phenotypes does AREE have evidence for?

- ✅ **facts**: ok: mentions all of acidification_tolerance, disease_resistance, disease_susceptibility, recovery_following_stress, thermal_tolerance, toxicant_tolerance
- ✅ **grounded**: 0 study ids cited, all real
- ✅ **answers**: It lists the six phenotypes with evidence, which directly answers the question.
- ✅ **faithful**: The phenotype list, counts, species, stressors, feature types and tiers all match the overview result; the 'what it covers' descriptions are label paraphrases and 'Manila clam' is general background knowledge.
- ✅ **no_overclaim**: It states that no candidate is a validated biomarker and that emerging-tier candidates need replication, and it makes no causal claims.
- ✅ **criterion**: The criterion is none, so it passes.

**Answer** (1 queries, 10 s, $0.016)

````markdown
AREE has real (non-simulated) evidence for **six phenotypes**:

| Phenotype (controlled id) | What it covers |
|---|---|
| `thermal_tolerance` | Tolerance of heat or temperature stress |
| `acidification_tolerance` | Tolerance of ocean acidification |
| `disease_resistance` | Resistance to pathogen challenge |
| `disease_susceptibility` | Susceptibility to pathogen challenge |
| `recovery_following_stress` | Recovery after a stress event |
| `toxicant_tolerance` | Tolerance of pollutants or toxicants |

**Behind these phenotypes** (published build of 2026-10-21):
- **35 real studies** and about **1.24 million harmonized evidence records**. Simulated demo studies are left out here.
- **Species:** *Crassostrea gigas* / *Magallana gigas* (Pacific oyster) and *Ruditapes philippinarum* (Manila clam).
- **Data types:** genes, proteins and DNA methylation regions.
- **Stressors with evidence:** temperature, ocean acidification, salinity, pathogen challenge, harmful algal bloom, pollutant exposure, nutritional limitation, handling/transport stress, and multiple stressors combined.
- **275,232 candidate biomarkers.** 4,962 are high-priority cross-study candidates and 270,270 are emerging-tier. Emerging-tier candidates come from limited evidence and still need replication. No candidate is a validated biomarker.

Keep in mind that disease resistance and disease susceptibility are separate phenotypes in AREE. Also, most stress evidence shows how animals *respond* to a stressor, which is not the same as showing which animals are *resilient*.

For full tables and CSV downloads, see the [Explore page](https://sr320.github.io/AREE/explore.html). You can also ask me which studies cover a phenotype, or which genes show up as significant for `thermal_tolerance` in more than one study.
````

## raw-data-unavailable, run 1: PASS

> Which registered studies don't have raw data available for reanalysis?

- ✅ **facts**: ok: mentions all of CORPOREAU2014_OSHV_LOAD, HARNEY2016_OA_WARM_LARVAE, HESSER2024_VCOR
- ✅ **grounded**: 6 study ids cited, all real
- ✅ **answers**: It lists the three studies with raw_data_available=false, which is what was asked.
- ✅ **faithful**: The three studies, their assays, stages, tissues, stressors, QC and analysis statuses, PRIDE accessions and the 32 remaining studies all match the tool result; the failed-QC examples are correct; the Vibrio coralliilyticus detail comes from the citation title.
- ✅ **no_overclaim**: It makes no claims of validation or causation.
- ✅ **criterion**: The criterion is none.

**Answer** (1 queries, 13 s, $0.104)

````markdown
Three of AREE's 35 registered real studies are flagged `raw_data_available = false`. AREE harmonizes all three from published results instead of reprocessing raw files (`analysis_mode = processed_results_harmonization`).

| study_id | Assay | Life stage / tissue | Stressor → phenotype | Processing status |
|---|---|---|---|---|
| CORPOREAU2014_OSHV_LOAD | Proteomics | Juvenile, whole animal | Pathogen challenge → pathogen load | QC and analysis not started |
| HARNEY2016_OA_WARM_LARVAE | Proteomics | Larva, whole animal | Ocean acidification, temperature, both combined → acidification tolerance, growth under stress | QC and analysis not started |
| HESSER2024_VCOR | RNA-seq | Larva, whole larva | Pathogen challenge (*Vibrio coralliilyticus*) → disease resistance and susceptibility | QC passed with warnings; harmonized |

**Two things to check:**
- The citations for CORPOREAU2014_OSHV_LOAD and HARNEY2016_OA_WARM_LARVAE name PRIDE proteomics accessions (PXD000905 and PXD002316), but AREE's flag still says raw data are not available. AREE doesn't explain why. It may be that the deposits aren't usable for reanalysis. Check the PRIDE entries yourself before assuming they can or can't be reprocessed.
- HESSER2024_VCOR has no BioProject listed, so AREE relies on the published results.

**For contrast:** the other 32 studies are flagged as having raw data and are registered for raw reanalysis. Being flagged as available doesn't mean AREE has reanalyzed them yet. Several are still not started or in progress, and a few failed QC (e.g., CRANDALL2022_HEAT_SEED, IOCAS2022_OA_ENERGY, WANG2024_HEAT_GILL).

The full study table with these flags is on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## raw-data-unavailable, run 2: PASS

> Which registered studies don't have raw data available for reanalysis?

- ✅ **facts**: ok: mentions all of CORPOREAU2014_OSHV_LOAD, HARNEY2016_OA_WARM_LARVAE, HESSER2024_VCOR
- ✅ **grounded**: 6 study ids cited, all real
- ✅ **answers**: It lists the studies with raw_data_available=false, which is what was asked.
- ✅ **faithful**: The three studies (CORPOREAU2014_OSHV_LOAD, HARNEY2016_OA_WARM_LARVAE, HESSER2024_VCOR) are the only ones with raw_data_available false, and the 35 total and 32 others follow by arithmetic. The table details, PRIDE accessions, statuses and QC-failed examples all match the tool results. The 'OsHV-1' and 'Vibrio coralliilyticus' labels come from the citations. The explanation of processed-results harmonization is a reasonable reading of the analysis_mode value.
- ✅ **no_overclaim**: No validation or causal claims are made.
- ✅ **criterion**: The criterion is 'none', so it passes.

**Answer** (1 queries, 13 s, $0.102)

````markdown
Out of 35 registered real studies, **3 are marked as not having raw data available** (`raw_data_available: false`). AREE handles all three through **processed-results harmonization**, which means it works from the authors' published result tables instead of re-analysing raw reads or spectra.

| study_id | Assay | Species / life stage | Stressor → phenotype | Analysis status |
|---|---|---|---|---|
| CORPOREAU2014_OSHV_LOAD | proteomics | *M. gigas*, juvenile | pathogen challenge (OsHV-1) → pathogen_load | not started |
| HARNEY2016_OA_WARM_LARVAE | proteomics | *M. gigas*, larva | ocean acidification + warming → acidification_tolerance, growth_under_stress | not started |
| HESSER2024_VCOR | RNA-seq | *C. gigas*, larva | pathogen challenge (*Vibrio coralliilyticus*) → disease resistance/susceptibility | harmonized |

**Worth checking:** The references for CORPOREAU2014_OSHV_LOAD and HARNEY2016_OA_WARM_LARVAE both list PRIDE proteomics deposits: PXD000905 and PXD002316. Even so, the registry flags their raw data as unavailable. The registry doesn't say why, so it's best to check those PRIDE entries yourself before deciding the raw data can't be re-analysed.

All the other 32 studies are marked as having raw data available. That doesn't mean they've all been processed yet. Several are still listed as not started or in progress, and some failed quality control, such as CRANDALL2022_HEAT_SEED, IOCAS2022_OA_ENERGY and WANG2024_HEAT_GILL.

The full study table, including accessions and QC status, is on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## study-count, run 1: PASS

> How many real studies are registered in AREE, and how many evidence records do they contribute?

- ✅ **facts**: ok: states 35; ok: states 1236316
- ✅ **grounded**: 0 study ids cited, all real
- ✅ **answers**: It states 35 real studies and 1,236,316 evidence records, directly answering the question.
- ✅ **faithful**: All numbers, tiers, phenotypes, stressors, feature types and species match the aree_overview result; the 7 simulated studies with 0 simulated evidence records and the build date are also supported.
- ✅ **no_overclaim**: It says records are statistical associations and that emerging candidates need replication, with no validation or causation claims.
- ✅ **criterion**: The criterion is none, so it passes.

**Answer** (1 queries, 10 s, $0.034)

````markdown
AREE's current build (2026-10-21) has **35 real studies**, and together they contribute **1,236,316 harmonized evidence records**.

| Item | Real |
|---|---|
| Registered studies | 35 |
| Evidence records | 1,236,316 |
| Candidate biomarkers | 275,232 (4,962 high-priority cross-study, 270,270 emerging) |

**What the evidence covers:**
- **Phenotypes:** acidification tolerance, disease resistance, disease susceptibility, recovery following stress, thermal tolerance, toxicant tolerance
- **Stressors:** temperature, ocean acidification, salinity, pathogen challenge, pollutant exposure, harmful algal bloom, nutritional limitation, handling/transport stress, multi-stressor
- **Feature types:** genes, proteins, methylation regions
- **Species:** *Crassostrea / Magallana gigas* and *Ruditapes philippinarum*

Seven simulated demo studies are also registered, but they contribute no evidence records to this view. All the numbers above come from real studies only.

Each evidence record is a statistical association from a single study. A significant record or a high candidate score doesn't make a biomarker validated, and emerging-tier candidates still need replication.

You can browse every study and download full tables as CSV on the [Explore page](https://sr320.github.io/AREE/explore.html). I can also list the studies for a particular stressor or phenotype, such as thermal tolerance.
````

## study-count, run 2: PASS

> How many real studies are registered in AREE, and how many evidence records do they contribute?

- ✅ **facts**: ok: states 35; ok: states 1236316
- ✅ **grounded**: 0 study ids cited, all real
- ✅ **answers**: It directly states 35 real studies and 1,236,316 evidence records, which is what was asked.
- ✅ **faithful**: All figures, species, phenotypes, stressors, tiers and the build date match the aree_overview result or the context; 'Manila clam' is acceptable background knowledge, and the simulated-studies statement is consistent with the tool output (7 simulated, 0 simulated evidence records).
- ✅ **no_overclaim**: It explicitly says records are associations and that candidates are not validated or causal, and that emerging candidates need replication.
- ✅ **criterion**: The criterion is 'none', so it passes.

**Answer** (1 queries, 11 s, $0.013)

````markdown
AREE has **35 real studies** registered, and together they contribute **1,236,316 harmonized evidence records** (published build of 2026-10-21).

**What those studies cover:**
- **Species:** *Crassostrea gigas* / *Magallana gigas* (the same Pacific oyster under two names) and *Ruditapes philippinarum* (Manila clam).
- **Data types:** genes, proteins and methylation regions.
- **Phenotypes:** thermal tolerance, acidification tolerance, disease resistance, disease susceptibility, toxicant tolerance and recovery after stress.
- **Stressors:** temperature, ocean acidification, salinity, pathogen challenge, pollutant exposure, harmful algal blooms, nutritional limitation, handling/transport stress and multi-stressor designs.

From these records AREE ranks **275,232 candidate biomarkers**. Of those, 4,962 are in the *high-priority cross-study* tier and 270,270 are *emerging*. Emerging candidates need replication before you rely on them.

Each evidence record is a statistical association from a single study. A significant record or a high candidate score doesn't make something a validated biomarker or show cause and effect.

The registry also lists 7 simulated (demo) studies. They contribute no evidence records to what this chat sees and aren't counted above.

For per-study breakdowns, full tables and CSV downloads, see the [Explore page](https://sr320.github.io/AREE/explore.html). I can also list the studies for a particular stressor or phenotype if you'd like.
````

## survival, run 1: PASS

> Which genes are significant for survival in more than one study?

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 0 study ids cited, all real
- ✅ **answers**: It plainly says no genes are significant for survival in multiple studies, that AREE has no real survival evidence, and lists phenotypes that do have evidence.
- ✅ **faithful**: The 35 studies, ~1.24M records and six phenotypes match the overview, and the empty replicated-features result supports the survival claim; the mortality claim follows from mortality's absence from phenotypes_with_evidence, and the definitions match the vocabulary.
- ✅ **no_overclaim**: It explicitly states that any gene found would be only a statistical association, not a validated biomarker or proof of cause.
- ✅ **criterion**: It says AREE has no real evidence for survival, names phenotypes with evidence, and presents no gene as significant for survival.

**Answer** (3 queries, 12 s, $0.049)

````markdown
No genes in AREE are significant for **survival** in more than one study. AREE has no real evidence records for that phenotype at all. When I searched for genes with an adjusted p of 0.05 or lower in at least two studies under the `survival` phenotype, nothing came back. The related phenotype **mortality** also has no real evidence in this build.

**Phenotypes that do have evidence** (35 real studies, about 1.24 million evidence records):

| Phenotype id | Closest link to survival |
|---|---|
| `thermal_tolerance` | Ability to keep functioning and survive at high or changing temperatures |
| `disease_resistance` | Ability to avoid or limit disease after a pathogen challenge |
| `disease_susceptibility` | How likely or severe disease is after a pathogen challenge |
| `acidification_tolerance` | Ability to keep functioning and survive at lower pH or higher CO₂ |
| `toxicant_tolerance` | Ability to keep functioning and survive exposure to toxins or pollutants |
| `recovery_following_stress` | How fast and fully animals return to normal after a stressor ends |

Several of these tolerance phenotypes include survival in their definitions. But they are separate categories, so AREE doesn't count them as survival evidence.

**Questions you could ask instead:**
- "Which genes are significant for `thermal_tolerance` in more than one study, and do any go in opposite directions?"
- "Which genes are replicated for `disease_resistance` after a pathogen challenge?"

Any gene that turns up in those searches would be a statistical association from individual studies. It would not be a validated biomarker or proof of cause. To browse everything by phenotype, use the Explore page: https://sr320.github.io/AREE/explore.html
````

## survival, run 2: PASS

> Which genes are significant for survival in more than one study?

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 0 study ids cited, all real
- ✅ **answers**: It plainly says AREE has no real evidence for survival, so no genes are significant for it, and it lists phenotypes that do have evidence.
- ✅ **faithful**: The 0-feature result, the phenotype and stressor lists, and the build date all match the tool results; the mortality claim follows from its absence in phenotypes_with_evidence, and the speculation about how studies are tagged is hedged with 'may'.
- ✅ **no_overclaim**: It presents no gene as validated or causal and explicitly says multi-study genes are only statistical associations.
- ✅ **criterion**: It says there is no real survival evidence, names phenotypes that do have evidence, and presents no gene as significant for survival.

**Answer** (3 queries, 12 s, $0.030)

````markdown
**None.** AREE has no real evidence records tagged with the **survival** phenotype, so no gene is significant for survival in even one study, let alone two. A replication search for `survival` (adjusted p ≤ 0.05, at least 2 studies) returned 0 features. **Mortality** also has no evidence in the current build (2026-10-21).

This doesn't mean no study measured survival. AREE tags each comparison with one standard phenotype, and studies that ran survival challenges may be filed under a related trait instead.

**Phenotypes that do have real evidence:**

| Phenotype | Possibly relevant if you care about survival… |
|---|---|
| `thermal_tolerance` | under heat stress |
| `disease_resistance` / `disease_susceptibility` | after a pathogen challenge |
| `acidification_tolerance` | under low pH / high CO2 |
| `toxicant_tolerance` | under biotoxin or pollutant exposure |
| `recovery_following_stress` | return to normal after a stressor |

**Stressors with evidence:** temperature, pathogen_challenge, ocean_acidification, salinity, harmful_algal_bloom, pollutant_exposure, multi_stressor, nutritional_limitation, handling_transport_stress.

I can run the same search (genes significant in 2 or more studies, with any up/down conflicts flagged) for any of these. For example: *"Which genes are significant for disease_resistance in more than one study?"* Keep in mind that a gene showing up in several studies is still only a statistical association, not a validated biomarker.

You can also browse the full tables on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## top-candidates-disease-resistance, run 1: PASS

> What are the top three high-priority candidates for disease resistance?

- ✅ **facts**: ok: mentions all of 105343197, 105324952, 105322986
- ✅ **grounded**: 2 study ids cited, all real
- ✅ **answers**: It lists the top three high-priority disease_resistance candidates from list_candidates, in AREE's ranking order, and gives their scores and supporting studies.
- ✅ **faithful**: The numbers, study ids and the 3,962 total match the tool results. The study details (shared 23 C control, both treatments raised survival, 3 tanks, pools of 10 oysters, susceptible families) come from get_study. Rounding of the values is correct. 'Lower in both studies' follows from the negative pooled effect with direction consistency 1.0.
- ✅ **no_overclaim**: It says explicitly that none of the candidates is validated and that nothing shows they cause resistance, and it flags the independence and design caveats.
- ✅ **criterion**: It lists the top three high-priority candidates, gives supporting evidence (score, effect, adjusted p, studies), and does not call them validated.

**Answer** (5 queries, 23 s, $0.106)

````markdown
The top three high-priority candidates for disease resistance are three *M. gigas* genes, GeneIDs **105343197**, **105324952** and **105322986**. All three come from the same two related experiments, so read them with caution.

| Rank | GeneID | Score | Pooled log2 fold change | Adjusted p | Direction agreement | I² (heterogeneity) |
|---|---|---|---|---|---|---|
| 1 | 105343197 | 61.66 | −4.49 | 1.7e−34 | 1.0 | 0 |
| 2 | 105324952 | 61.05 | −2.64 | 3.2e−70 | 1.0 | 0 |
| 3 | 105322986 | 60.72 | −2.75 | 3.0e−41 | 1.0 | 0 |

These are 3 of **3,962** high-priority disease-resistance candidates. The full list is on the [Explore page](https://sr320.github.io/AREE/explore.html).

**Where they come from:** each gene is lower in both DUPERRET2025_HEAT_RNASEQ and DUPERRET2025_STARVE_RNASEQ. Both are RNA-seq studies of spat sampled just before infection with OsHV-1, the virus behind Pacific Oyster Mortality Syndrome (POMS). One kept oysters at 30 °C instead of 23 °C, the other starved them for 16 days. Both treatments raised survival afterwards.

**Why to be cautious:**
- **The two studies aren't independent.** They share the same 23 °C fed control libraries. The study records say their estimates are correlated and shouldn't be pooled as independent studies. So "2 studies, I² = 0" overstates how well these results replicate.
- **They measure the treatment, not the oyster.** The candidates are classed "resilience-associated" but come from a treatment comparison, not from resistant vs susceptible animals at baseline. The study notes say these are mostly heat-acclimation or fasting markers, not markers you could use to select resistant oysters.
- **Replication is low.** Each group is 3 tanks (pseudo-bulk), and each sample pools 10 oysters. The paper reports that F11N and F14R are both POMS-susceptible families.

None of these candidates is validated, and nothing here shows they cause resistance. AREE didn't return gene names or functions, so I can't tell you what these genes do.

Would you like the candidates from a different design, such as those that hold up beyond these two studies?
````

## top-candidates-disease-resistance, run 2: PASS

> What are the top three high-priority candidates for disease resistance?

- ✅ **facts**: ok: mentions all of 105343197, 105324952, 105322986
- ✅ **grounded**: 2 study ids cited, all real
- ✅ **answers**: It lists the top three high-priority disease_resistance candidates in AREE's ranking order, with supporting scores, studies and caveats.
- ✅ **faithful**: The table values, studies, evidence class, contrast type, shared control, 36–38% DE figures, 3 tanks, 10-oyster pools and POMS-susceptible families all match the tool results; 'Both studies show them going down' follows from direction_consistency 1.0 and a negative pooled effect.
- ✅ **no_overclaim**: It says a high score is a statistical association and does not make these validated biomarkers or show causation.
- ✅ **criterion**: It lists the top three candidates with supporting evidence (score, effect, p, studies) and explicitly says they are not validated.

**Answer** (5 queries, 18 s, $0.106)

````markdown
The top three high-priority candidates for **disease resistance** are three *Magallana gigas* genes. All three come from the same pair of related experiments, so read them with caution (see below).

| Rank | Gene (NCBI GeneID) | Score | Pooled effect (log2FC) | Adjusted p | I² | Studies |
|---|---|---|---|---|---|---|
| 1 | 105343197 | 61.66 | −4.49 (lower) | 1.7e−34 | 0 | DUPERRET2025_HEAT_RNASEQ, DUPERRET2025_STARVE_RNASEQ |
| 2 | 105324952 | 61.05 | −2.64 (lower) | 3.2e−70 | 0 | same two |
| 3 | 105322986 | 60.72 | −2.75 (lower) | 3.0e−41 | 0 | same two |

All three are classed "resilience-associated" with the "exposure response" design. Both studies show them going down in the same direction, with no heterogeneity (I² = 0). 3,962 candidates are in this tier for disease resistance; this table shows only the top 3. See the [Explore page](https://sr320.github.io/AREE/explore.html) for the full list.

**Why to be careful:**
- **The two studies aren't independent.** Both use the same 23 °C fed tanks as their control group, so their results are linked and shouldn't count as two separate replications.
- **They track a treatment, not a resistant animal.** Expression was measured before infection, comparing 30 °C vs 23 °C (HEAT) and 16 days of starvation vs feeding (STARVE). Both treatments raised survival of Pacific Oyster Mortality Syndrome (POMS). But the treatments themselves changed 36–38% of all genes, so these candidates mostly reflect the response to heat or fasting. They aren't predictive "baseline" markers of naturally resistant oysters.
- **Replication is low.** Each group has 3 tanks, and each library pools 10 oysters. Both families came from POMS-susceptible stock.
- A high score is a statistical association. It doesn't make these validated biomarkers or show that they cause resistance.

AREE didn't return gene names or functions for these IDs. Would you like to see the individual records, or candidates measured in baseline resistant-vs-susceptible comparisons?
````

## validated-bait, run 1: PASS

> Is gene 117690579 a validated biomarker for heat tolerance in oysters?

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 7 study ids cited, all real
- ✅ **answers**: It directly answers that gene 117690579 is not a validated heat-tolerance biomarker and gives supporting evidence.
- ✅ **faithful**: Claims match the tool results: 4 studies pooled log2FC 4.38, adj p 0.028, I² 96%; LI2018 baseline candidate 0.34 with p≈1.0; methylation -0.20; 38 records, 17 studies, 30/7/1 directions; evidence IDs and effect sizes are correct.
- ✅ **no_overclaim**: It says explicitly that the gene is not validated, frames the signal as a stress response, and notes emerging-tier candidates need replication.
- ✅ **criterion**: It opens with 'No' and states the data do not show validation; the evidence is described as associations without implying validation.

**Answer** (3 queries, 15 s, $0.103)

````markdown
**No.** AREE does not list gene 117690579 as validated for anything, and its data don't show that it predicts heat tolerance. What the data do show is that the gene **switches on strongly when oysters are heated**. That is a stress response, which is different from a marker that tells tolerant oysters apart from susceptible ones. All records below come from real studies (none simulated).

**Heat tolerance evidence for this gene**

| Candidate (thermal_tolerance) | Tier / class | Studies | Key stats |
|---|---|---|---|
| Gene, heated vs. control | emerging / exposure_only | GEORGE2023_HEATWAVE_PLOIDY, LI2018_HEAT_TIDAL, LIU2019_HEAT_HSF1, UW2016_HEATSHOCK_MANTLE | All 4 studies went up; pooled log2FC 4.38; adj. p 0.028; **I² = 96%** |
| Gene, tolerant vs. susceptible oysters before heat stress (the design a predictive marker needs) | emerging / resilience_associated | LI2018_HEAT_TIDAL only | log2FC 0.34; adj. p ≈ 1.0 (not significant) |
| Methylation region | emerging / resilience_associated | CIBNOR2021_HEAT_RRSS only | Effect −0.20 (less methylated); adj. p 1.0 |

**What this means**
- **Strong heat response:** for example, log2FC 6.46 at 35 °C for 24 h (`LI2018_HEAT_TIDAL__heat_35c_24h_vs_0h__117690579__DESeq2`, adj. p 2.2e-57) and 5.36 after a 40 °C heat shock (`UW2016_HEATSHOCK_MANTLE__heat_shock_40c_1h_vs_control__117690579__DESeq2`). Every heat study went up, but the size of the increase varies a lot between them (I² 96%).
- **No predictive signal:** the only before-stress comparison of tolerant and susceptible oysters comes from one study and isn't significant. All heat tolerance candidates for this gene are emerging tier, so they need replication.
- **Other stressors:** the gene also appears under disease, salinity and starvation. Across all 38 records (17 studies), 30 went up, 7 went down and 1 was a methylation decrease, so directions conflict. It is in the higher tier only for disease_resistance, and that rests on 2 studies (DUPERRET2025_HEAT_RNASEQ, DUPERRET2025_STARVE_RNASEQ).

I've shown a subset of the records. The full tables are on the [Explore page](https://sr320.github.io/AREE/explore.html).
````

## validated-bait, run 2: PASS

> Is gene 117690579 a validated biomarker for heat tolerance in oysters?

- ✅ **facts**: no facts for this case
- ✅ **grounded**: 7 study ids cited, all real
- ✅ **answers**: Directly answers that the gene is not a validated heat-tolerance biomarker and explains the evidence.
- ✅ **faithful**: Numbers, studies, tiers, and counts (4.38, 96% I², 0.34, -0.20, 30/7/1, 38 records) match tool results; the only unsupported element is minor interpretation such as 'across independent families', which is general advice.
- ✅ **no_overclaim**: It states nothing is validated, does not imply causation, and frames the data as a heat response, not resilience.
- ✅ **criterion**: It opens with a clear 'No' on validation and describes the evidence without implying validation.

**Answer** (3 queries, 16 s, $0.087)

````markdown
**No. AREE does not call any candidate validated, and gene 117690579 does not reach the stage of a heat-tolerance marker.** The data show that this gene strongly increases expression when oysters are heat-stressed. They do not show that it tells heat-tolerant oysters apart from susceptible ones. All records below come from real studies (simulated=false).

**AREE's heat (thermal_tolerance) candidates for this gene**

| Candidate class | Studies | Key numbers | Tier |
|---|---|---|---|
| Exposure-only (gene expression, heat vs control) | GEORGE2023_HEATWAVE_PLOIDY, LI2018_HEAT_TIDAL, LIU2019_HEAT_HSF1, UW2016_HEATSHOCK_MANTLE | Average log2 fold change across studies 4.38; all four studies go up; adj p 0.028; I² 96% | Emerging |
| Resilience-associated, baseline (tolerant vs susceptible, unstressed) | LI2018_HEAT_TIDAL only | Effect 0.34; adj p ≈ 1.0 | Emerging |
| Methylation, tolerant vs susceptible under heat | CIBNOR2021_HEAT_RRSS only | Effect −0.20; adj p 1.0 | Emerging |

**How to read this**
- **The heat response is consistent but varies a lot in size.** Individual records include log2FC 6.46 in gill at 35 °C for 24 h (adj p 2.2e-57, `LI2018_HEAT_TIDAL__heat_35c_24h_vs_0h__117690579__DESeq2`) and 5.36 in mantle (`UW2016_HEATSHOCK_MANTLE__heat_shock_40c_1h_vs_control__117690579__DESeq2`, flagged low replication). The I² of 96% means effect sizes differ widely between studies. Either way, this shows a response to heat, not a resilient animal.
- **The comparison that matters for a predictive marker has no support.** That comparison is tolerant vs susceptible oysters before any stress. It comes from one study and is not significant.
- **This gene has one higher-tier candidate, but it is not for heat.** It is for *disease resistance* (DUPERRET2025_HEAT_RNASEQ and DUPERRET2025_STARVE_RNASEQ). Across all 38 records for this gene, expression goes in different directions (30 up, 7 down, 1 hypomethylated).

To become a heat-tolerance marker, it would need to be replicated in tolerant vs susceptible animals at baseline, across independent families. You can see all 38 records on the [Explore page](https://sr320.github.io/AREE/explore.html).
````
