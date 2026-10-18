# Contrast types: baseline differences vs. responses

`resilience_classification` says what a comparison measured: a resilience
outcome, a disease response, a stress response, or exposure alone
([resilience_vs_exposure.md](resilience_vs_exposure.md)). It does not say
*between whom* or *when*. Two comparisons can both be resilience evidence and
still ask different questions:

- Do intertidal oysters, which tolerate heat better, express a gene differently
  from subtidal oysters **before either is heated**?
- How does that gene respond **when an oyster is heated**?

A breeder can only use the first kind directly. A marker that already differs
between resilient and susceptible animals under ordinary conditions can be read
before any stress, on animals that are then kept. A gene that responds to heat
tells you an animal is being heated.

So every comparison also declares a `contrast_type`
(`comparisons[].contrast_type` in `schemas/study.schema.json`, required).

## The four contrast types

| `contrast_type` | The two arms are | Sampled | Example |
|---|---|---|---|
| `baseline_phenotype_difference` | groups that differ in resilience **by origin**: line, family, population, ploidy or age | under common, unstressed conditions | LI2018 intertidal vs. subtidal populations at 0 h |
| `phenotype_difference_under_exposure` | groups that differ in resilience, by origin, outcome or prior conditioning | during or after a challenge | CIBNOR2021 heat-resistant vs. heat-susceptible families during heat; CORPOREAU2014 high vs. low viral load |
| `exposure_response` | a treatment and its control | any time | heat vs. control; infected vs. uninfected; starved vs. fed |
| `phenotype_by_exposure_interaction` | the difference between two groups' responses to one exposure | in a 2 × 2 design | resistant line's heat response minus susceptible line's heat response |

How to choose:

1. **Do the two arms differ by something done to them?** If so, it is an
   `exposure_response`, even when sampling happens before a later challenge.
   DUPERRET2025's starved vs. fed oysters, sampled before infection, are an
   exposure response: starvation made them more resistant, but the contrast is
   a diet treatment, not a difference between resilient and susceptible animals.
2. **Do the arms differ in resilience by origin, and were both sampled under
   the same unstressed conditions?** That is a `baseline_phenotype_difference`.
   The resilience difference must be measured, in the study or in a cited
   source, not assumed; record where in `limitations`. A baseline contrast is
   resilience evidence by construction, so validation warns if its
   `resilience_classification` is anything but `resilience`.
3. **Do the arms differ in resilience but were sampled during or after the
   challenge?** That is a `phenotype_difference_under_exposure`. A difference
   here can be a consequence of the outcome rather than a cause: dying animals
   express different genes.

When in doubt, choose `exposure_response`. It is the type that claims least.

## What the pipeline does with it

- **Pooling.** `contrast_type` is part of the meta-analysis grouping key and of
  each Benjamini-Hochberg test family. A baseline contrast never pools with an
  exposure response, even under the same phenotype and feature. Cross-species
  pooling through orthogroups keeps the same split.
- **Primary comparisons.** `meta_analysis_primary` is chosen within each
  contrast type. A study whose primary comparison is a heat response still
  contributes its baseline contrast to the baseline pool. Before contrast types,
  LI2018's intertidal-vs-subtidal contrast was set aside because the study's
  primary comparison was the 35 °C response.
- **Multi-omics convergence** counts layers significant under the same
  phenotype **and** contrast type. A methylation difference at baseline and an
  expression response to heat are two questions, not convergence.
- **Ranking.** Every candidate carries its `contrast_type`. A candidate whose
  contrast type is `baseline_phenotype_difference` and whose evidence class is
  `resilience_associated` is flagged `is_predictive_baseline`. The flag does not
  change its tier or score: tier still says how well replicated it is.
- **Reporting.** `aree top-candidates` opens with a "Predictive baseline
  candidates" section. Evidence cards state the contrast type and what it does
  and does not show. A card for a non-default contrast type has the type as a
  filename suffix, so it never overwrites the exposure-response card for the
  same feature, and its per-study table names each comparison, so the sign of
  the effect can be read against the arms. The query layer and MCP server can
  filter candidates by `contrast_type`.
- **Website.** The [home page](index.qmd) has a "Predictive baseline
  candidates" section: the significant ones, by tier then score, with gene
  names and links to their evidence cards, and a plain statement when all rest
  on a single study. On the [Explore page](explore.qmd), the Candidates tab
  filters by contrast type and has a "Predictive baseline only" toggle; the
  candidates table carries `gene_symbol` and `gene_description`, and
  `is_predictive_baseline` is a true/false column usable directly in SQL.
- **Old evidence.** An evidence table harmonized before this field existed has
  no `contrast_type`. The meta-analysis refuses it and names the studies to
  re-harmonize, rather than defaulting them to `exposure_response`.

## What replication means for a baseline candidate

Groups that differ in origin differ in more than resilience: genetic
background, growth rate, the site they came from. One resistant-vs-susceptible
pair cannot tell a resilience marker from a background difference. Treat
`k_studies` for a baseline candidate as the number of independent group pairs
that agree, and look for:

1. the same difference in **independent** line, family or population pairs
   (the high-priority tier's two-study gate);
2. ideally, **individual-level** validation: sample animals non-lethally before
   a challenge, record who survives, and test whether the marker predicts it.

The second is the step that turns a predictive candidate into a marker.

## Registered baseline contrasts

| Study | Comparison | Status |
|---|---|---|
| LI2018_HEAT_TIDAL | `intertidal_vs_subtidal_0h` | harmonized |
| DECICOMP2025_POMS_AGE | `age28mo_vs_age4mo_basal` (primary), `age16mo_vs_age4mo_basal` | not started |
| DECICOMP2025_POMS_AGE_EMSEQ | `age28mo_vs_age4mo_basal` (primary), `age16mo_vs_age4mo_basal` | not started |

The DECICOMP2025 contrasts compare ages, and age differs in development as well
as POMS resistance; their `limitations` should say so when they are run.
