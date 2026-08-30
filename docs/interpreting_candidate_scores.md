# Interpreting candidate biomarker scores

AREE assigns every candidate (one row per pooled feature/phenotype/feature-type
group) a numeric score and a tier. This page explains both — and, more
importantly, explains that the tier is decided by hard gates checked
**before** the score, not by the score itself. A high score cannot promote a
weakly-supported candidate into a strong tier.

Candidate identity also includes `simulated` and `species_taxid`. Assay
diversity, mapping confidence, quality flags, forest data, and card filenames
are calculated within that partition, so a demo record or another species can
never promote or overwrite a real single-species candidate.

## The score: a transparent weighted sum

`src/prioritize/scoring.py` computes a 0-1 value for each named component and
combines them as a plain weighted sum, minus a heterogeneity penalty — there
is no hidden model or machine-learning step. Weights sum to 100 across the
positive components:

| component | weight | what it measures | source |
|---|---|---|---|
| `n_studies_score` | 20 | independent replication (`k_studies / 5`, capped at 1) | count of distinct `study_id` |
| `direction_consistency_score` | 20 | agreement in effect direction across studies | meta-analysis `direction_consistency` |
| `sample_size_score` | 10 | total biological replication (`total_sample_size / 100`, capped) | sum of `sample_size` |
| `effect_magnitude_score` | 10 | typical absolute pooled effect (`|pooled_effect| / 2`, capped) | pooled `effect_size` |
| `significance_score` | 10 | adjusted-significance strength (`-log10(p) / 5`, capped) | pooled `p_value` |
| `phenotype_relevance_score` | 10 | resilience (1.0) > disease (0.8) > stress_response (0.4) > exposure_only (0.1) | phenotype's `resilience_relevance` |
| `assay_diversity_score` | 10 | number of distinct molecular layers/feature types (`n_distinct_assays / 3`, capped) | distinct `feature_type` |
| `context_breadth_score` | 5 | spread across tissues × life stages (`/3`, capped) | distinct `tissue` × `life_stage` |
| `mapping_confidence_score` | 3 | trust in identifier harmonization (worst mapping among contributing records) | `mapping_confidence` |
| `quality_score` | 2 | study/data quality (`1 - n_quality_flags/5`, floored at 0) | `quality_flags` |

Separately, a **heterogeneity penalty** of up to 15 points is subtracted,
scaled by `i_squared / 100` — high cross-study heterogeneity directly reduces
the final score rather than being ignored.

```
score = sum(weight[c] * component[c] for c in components) - 15 * (i_squared / 100)
```

clipped to `[0, 100]`. The function (`candidate_score` in `scoring.py`) is
pure — identical inputs always produce the identical score — which is what
makes candidate-score reproducibility directly testable.

Notice `mapping_confidence_score` and `quality_score` are weighted quite
lightly (3 and 2 points respectively) relative to `n_studies_score` and
`direction_consistency_score` (20 each). The score alone will not sink a
candidate just for imperfect identifier mapping or a couple of quality
flags — which is exactly why those two properties are instead enforced as
**hard gates**, not left to the weighted sum, for the top tier (see below).

## The tiers, and the gates that decide them

`src/prioritize/rank.py` computes tier membership from hard boolean
conditions, evaluated independently of the numeric score:

```python
is_high_priority = (
    k_studies >= 2
    and phenotype_relevance_score > 0.1
    and direction_consistency >= 0.7
    and quality_score >= 0.4
)
is_multi_omics = n_distinct_assays >= 2
```

Tier assignment, in order:

1. **`high_priority_cross_study`** — requires **all** of:
   - `k_studies >= 2` (at least two independent studies)
   - `phenotype_relevance_score > 0.1` (not `exposure_only` — see
     [resilience_vs_exposure.md](resilience_vs_exposure.md))
   - `direction_consistency >= 0.7`
   - `quality_score >= 0.4` (i.e. fewer than ~3 of 5 quality flags present)
2. **`multi_omics_convergence`** — evidence from `>= 2` distinct
   `feature_type`s mapped to the same standardized feature (checked only if
   the candidate did not already qualify for `high_priority_cross_study`).
3. **`emerging`** — everything else: single-study support, or support that
   fails one or more of the high-priority gates. Always labeled as requiring
   replication.

**A high score cannot override these gates.** A candidate with `k_studies = 1`
that happens to have a huge effect size, tiny p-value, and perfect direction
consistency (trivially 1.0 with only one study) can still score numerically
high on the weighted sum — but it is placed in `emerging`, not
`high_priority_cross_study`, because it fails the `k_studies >= 2` gate.
Likewise, the conflicting-evidence example in
[interpreting_meta_analysis.md](interpreting_meta_analysis.md)
(`LOC105331241`/`sod1`, `direction_consistency = 0.5`) fails the
`direction_consistency >= 0.7` gate regardless of its score and cannot reach
`high_priority_cross_study`.

This is the concrete mechanism (see `design.md` section 6) that prevents
"statistically significant in one study" from reading as "validated": the
score is informative for ranking *within* a tier, but tier membership itself
is gated on replication, directional agreement, phenotype relevance, and
quality — properties a single striking p-value cannot substitute for.

## Reading the output

`aree build-evidence-cards` writes one evidence card per candidate under
`reports/evidence_cards/`, plus an `index.json` summarizing tier and score
for every candidate. Every candidate gets a card regardless of tier —
`emerging` candidates are not hidden, they are labeled. Each card shows the
component breakdown (`component_*` columns), the tier, the contributing
studies, and (where feasible) a forest-style summary of per-study effect
estimates, so a reviewer can see exactly which inputs drove the score rather
than trusting a single number.

## Related documentation

- [interpreting_meta_analysis.md](interpreting_meta_analysis.md)
- [resilience_vs_exposure.md](resilience_vs_exposure.md)
- [identifier_mapping.md](identifier_mapping.md)
- [design.md](design.md#6-candidate-scoring-framework-transparent-not-black-box)
