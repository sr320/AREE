# Interpreting meta-analysis output

`aree meta-analyze` pools effect sizes across studies for the same
standardized feature, phenotype, and feature type. This page explains each
output column and, critically, when a pooled number is a warning sign rather
than a summary you should trust.

## How pooling works

`src/meta_analysis/run.py` groups the evidence table by
`(feature_id_standardized, phenotype, feature_type)` — never across feature
types, since a log2FoldChange and a methylation percent-difference are not on
a comparable scale. Records with `mapping_confidence == "unresolved"` are
excluded from pooling (no stable identity to group on) but remain visible in
the raw evidence table. Within each group, effect sizes and standard errors
are pooled with a random-effects (DerSimonian–Laird) estimator implemented
directly in Python (`src/meta_analysis/pooling.py`) — see
[methods.md](methods.md) for the statistical detail.

## Output columns

| column | meaning |
|---|---|
| `k_studies` | number of distinct `study_id` values contributing |
| `studies` | pipe-separated list of contributing study IDs |
| `n_evidence_records` | number of evidence rows pooled (can exceed `k_studies` if a study contributes multiple comparisons) |
| `total_sample_size` | sum of `sample_size` across distinct (study, comparison) pairs |
| `pooled_effect` | random-effects pooled estimate of the effect size |
| `pooled_se` | standard error of the pooled effect |
| `ci_lower` / `ci_upper` | 95% CI on the pooled effect (`pooled_effect ± 1.96 × pooled_se`) |
| `z`, `p_value` | test statistic and two-sided p-value for the pooled effect |
| `q_statistic` | Cochran's Q, the weighted sum of squared deviations from the fixed-effect estimate |
| `i_squared` | percentage of total variation attributable to between-study heterogeneity rather than sampling error |
| `tau_squared` | estimated between-study variance component |
| `direction_consistency` | fraction of contributing records agreeing with the majority effect-direction sign |
| `distinct_tissues`, `distinct_life_stages` | breadth of biological contexts represented |
| `distinct_stressors` | pipe-separated list of standardized stressors contributing |
| `mapping_confidences` | union of mapping-confidence levels among contributing records |
| `quality_flags_union` | union of quality flags among contributing records |

## Reading `pooled_effect` and `ci_lower`/`ci_upper` correctly

A pooled effect with a CI that excludes zero, few studies, and low
heterogeneity is reasonable evidence of a consistent effect. A pooled effect
alone, without checking `i_squared` and `direction_consistency`, can be
actively misleading — averaging together a strong positive effect and a
strong negative effect can produce a pooled estimate near zero that looks
"unremarkable" instead of flagging a real conflict.

## Worked example: a genuine conflicting-evidence case

Running the demo meta-analysis for `larval_viability` genes
(`reports/meta_analysis/larval_viability_gene_meta_analysis.tsv`) produces
this row for gene `LOC105331241` (`sod1`, superoxide dismutase 1):

| field | value |
|---|---|
| `k_studies` | 2 |
| `studies` | `GIGAS_LARV05` \| `GIGAS_OA02` |
| `pooled_effect` | -0.191 |
| `pooled_se` | 0.800 |
| `ci_lower` / `ci_upper` | -1.759 / 1.377 |
| `p_value` | 0.811 |
| `i_squared` | 93.6% |
| `tau_squared` | 1.198 |
| `direction_consistency` | 0.5 |
| `mapping_confidences` | `exact` \| `inferred` |
| `quality_flags_union` | `identifier_mapping_uncertain` \| `low_replication` |

This is not "no effect." It is two studies reporting opposite-direction
effects for the same gene under the same phenotype (`larval_viability`):
`GIGAS_OA02` finds `sod1` (as `LOC105331241`) strongly **down**-regulated
under ocean acidification (log2FC ≈ -1.00, adjusted p = 0.03), while
`GIGAS_LARV05` finds it **up**-regulated under nutritional limitation
(log2FC ≈ 0.60). This is intentional in the demo data — the curation note in
`registry/studies/GIGAS_OA02.yaml` states the study "deliberately includes an
antioxidant-pathway gene (sod1) that is down-regulated here but up-regulated
in GIGAS_LARV05 under the same larval_viability phenotype, to demonstrate how
AREE surfaces conflicting same-phenotype cross-study evidence in
meta-analysis (low direction_consistency) rather than averaging it away
silently."

The signature to recognize:

- **`direction_consistency = 0.5`** with `k_studies = 2` means an exact
  50/50 split — the weakest possible directional agreement above "no
  majority."
- **`i_squared ≈ 93.6%`** means almost all of the variability across these
  two estimates is between-study heterogeneity, not sampling noise — the
  two studies are not measuring compatible effects.
- The wide CI (`-1.76` to `1.38`) spanning zero and both original
  study-level effect signs, combined with a non-significant pooled
  `p_value` (0.81), is the direct numerical consequence: pooling two
  opposite-direction effects returns "no discernible average effect,"
  which is a true statement about the pooled number but a misleading
  summary of the underlying evidence if read in isolation.

**How to read this correctly:** treat high `i_squared` (a common informal
threshold is >75%) together with `direction_consistency` well below 1.0 as a
signal that the two studies may be capturing genuinely different biology —
different stressors (`ocean_acidification` vs. `nutritional_limitation`),
different underlying mechanisms, or context-dependent regulation — not as
noise to be smoothed over by the pooled estimate. This is exactly why the
candidate-scoring framework does not let a low pooled p-value substitute for
directional agreement: `direction_consistency_score` is weighted as heavily
as `n_studies_score` (both 20 of 100 points), and heterogeneity is
subtracted as an explicit penalty — see
[interpreting_candidate_scores.md](interpreting_candidate_scores.md). This
candidate would score low on both counts and, separately, fails the
`high_priority_cross_study` hard gate (`direction_consistency >= 0.7`)
regardless of its overall score.

## Related documentation

- [interpreting_candidate_scores.md](interpreting_candidate_scores.md)
- [methods.md](methods.md)
- [design.md](design.md#6-candidate-scoring-framework-transparent-not-black-box)
