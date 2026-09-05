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
| `n_studies_score` | 20 | independent replication (`k / (k + 5)`) | count of distinct `study_id` |
| `direction_consistency_score` | 20 | agreement in effect direction across studies | meta-analysis `direction_consistency` |
| `sample_size_score` | 10 | total biological replication (`n / (n + 100)`) | sum of `sample_size` |
| `effect_magnitude_score` | 10 | absolute pooled effect (`\|e\| / (\|e\| + 2)`) | pooled `effect_size` |
| `significance_score` | 10 | adjusted-significance strength (`s / (s + 5)` where `s = -log10(q)`) | pooled `adjusted_p_value` |
| `phenotype_relevance_score` | 10 | resilience (1.0) > disease (0.8) > stress_response (0.4) > exposure_only (0.1) | phenotype's `resilience_relevance` |
| `assay_diversity_score` | 10 | molecular layers with a significant record for this feature under this phenotype (`n_supporting_layers / 3`, capped) | `feature_type` → `molecular_layer`, `adjusted_p_value` |
| `context_breadth_score` | 5 | spread across tissues × life stages (`c / (c + 3)`) | distinct `tissue` × `life_stage` |
| `mapping_confidence_score` | 3 | trust in identifier harmonization (worst mapping among contributing records) | `mapping_confidence` |
| `quality_score` | 2 | study/data quality (`1 - n_quality_flags/5`, floored at 0) | `quality_flags` |

### Why the unbounded components saturate instead of clipping

Every component whose input has no natural upper bound uses `v / (v + half)`,
which is strictly increasing everywhere, equals 0.5 at `half`, and approaches
but never reaches 1.0. Each `half` above is the value that used to be a hard
cap, so the ordering is unchanged and only the ceiling is gone.

This replaced `min(1, v / cap)` on **2026-09-05**, because clipping broke the
ranking exactly where it matters. On the first genome-wide real pool,
`|pooled_effect|` clipped at 2 and adjusted p clipped at `1e-5`; 126 of the
1,994 high-priority candidates therefore scored an identical **68.57** despite
pooled effects spanning 2 to 9.7 log2FC and adjusted p spanning `1e-5` to
`1e-75`. A reader choosing validation targets from the top of that list was
choosing between ties. After the change the largest tie group is 10, and those
are rounding collisions at two decimal places, not saturation.

The components with genuinely bounded inputs still clip: `direction_consistency`
and `quality_score` are 0-1 by construction, and `assay_diversity_score` is
capped by the number of molecular layers in the `feature_types` vocabulary.

Absolute scores are lower than before this change — nothing reaches a
component value of 1.0 any more. Compare candidates within one run, not across
this change.

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
    and adjusted_p_value <= 0.05
)
is_multi_omics = own_layer in supporting_layers and n_supporting_layers >= 2
```

Tier assignment, in order:

1. **`high_priority_cross_study`** — requires **all** of:
   - `k_studies >= 2` (at least two independent studies)
   - `phenotype_relevance_score > 0.1` (not `exposure_only` — see
     [resilience_vs_exposure.md](resilience_vs_exposure.md))
   - `direction_consistency >= 0.7`
   - `quality_score >= 0.4` (i.e. fewer than ~3 of 5 quality flags present)
   - `adjusted_p_value <= 0.05` — the pooled effect survives
     Benjamini–Hochberg control within its phenotype / feature-type family
     (see [interpreting_meta_analysis.md](interpreting_meta_analysis.md#multiple-testing)).
     Without this gate, two genome-wide studies promote every gene whose
     effects happen to share a sign: on the first real OsHV-1 pool that was
     12,658 of 23,094 genes, 8,533 of them with pooled `p > 0.05`.
2. **`multi_omics_convergence`** — at least two **molecular layers**
   (`transcriptomics`, `dna_methylation`, `proteomics`, `metabolomics`, as
   declared per feature type in
   `registry/controlled_vocabularies/feature_types.yaml`) each carry a
   significant record (study-level adjusted p ≤ 0.05) for the same
   standardized feature under the **same phenotype**, and the candidate's own
   layer is one of them. Checked only if the candidate did not already
   qualify for `high_priority_cross_study`.

   Three things this deliberately does *not* count: a DMR and a single-CpG
   DML for one gene (two views of one methylation layer); a gene expressed
   under heat and methylated under pathogen challenge (two questions, not
   convergent evidence for either); and two other layers converging on a gene
   the candidate itself shows no signal for. Before 2026-09-05 the gate was
   "≥ 2 distinct feature types anywhere for this gene", which promoted 27 of
   the 48 simulated demo candidates; under the current gate none of them
   qualifies, because every cross-layer overlap in the demo is
   cross-phenotype. The card still lists that adjacent evidence and marks
   each row with whether it counts.
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
is gated on replication, directional agreement, phenotype relevance,
quality, and family-wise significance — properties a single striking p-value
cannot substitute for. `significance_score` itself is also computed from
`adjusted_p_value`, not the raw pooled p.

## Evidence class: the tier answers a different question

The tier says **how well replicated** a signal is. It does not say **what the
signal is evidence of** — and a reader who sees `high_priority_cross_study`
will reasonably assume "resilience biomarker". For the first real pool that
assumption is wrong for every candidate.

Each candidate therefore carries two more columns, computed in
`src/prioritize/rank.py` and printed at the top of its evidence card:

| column | values | meaning |
|---|---|---|
| `evidence_class` | `resilience_associated`, `disease_associated`, `stress_response`, `exposure_only` | the `resilience_relevance` of the phenotype the contributing studies actually measured |
| `context_replication` | `single_study`, `single_context`, `multi_context` | whether the contributing studies vary in stressor, tissue, or life stage — replication of the *biology*, as distinct from replication of the study |

Both real OsHV-1 studies measure `disease_susceptibility`: animals were
challenged with the virus, but no survival or per-animal viral load was
recorded. So all 1,994 top-tier candidates are `disease_associated`, not
`resilience_associated`. They are well-replicated evidence of an infection
response. Whether any of them distinguishes an animal that survives from one
that does not is a question these data cannot answer, and the card says so
rather than leaving the tier to imply otherwise.

`context_replication` is `multi_context` for that pool, but only because the
two studies used different life stages (spat vs. juvenile); the stressor and
tissue are identical. The card names the dimensions that vary and the ones
that do not, instead of asserting a generality the pool has not established.

Because the phenotype gap outranks the replication gap, a card whose
`evidence_class` is not `resilience_associated` leads its recommended next
step with closing that gap — measuring the feature in animals whose outcome
was recorded — and only then gives the tier-appropriate validation advice.
Confirming an expression difference that was never tied to an outcome does not
turn it into a biomarker.

## Reading the output

`aree build-evidence-cards` ranks every candidate and writes two things under
`reports/evidence_cards/`:

- **`candidates.tsv`** — every ranked candidate, all tiers, with the full
  `component_*` breakdown, the meta-analysis fields, and a `card_file` column.
  This is the complete, auditable ranking; `emerging` candidates are not
  hidden, they are labeled.
- **One markdown card per candidate that carries a significant signal**, plus
  an `index.json` listing those cards. A candidate qualifies when its pooled
  `adjusted_p_value` is at or below `--max-adjusted-p` (default 0.05), *or*
  when at least one contributing study reported an adjusted p at or below it.
  The second clause is what keeps conflicting-direction candidates on a card:
  `LOC105331241`/`sod1` pools to `p = 0.81` because two studies disagree, but
  one of them reported `padj = 0.03`, and that conflict is exactly what a
  reader needs to see.

The threshold exists because a genome-wide pool yields one candidate per gene.
Rendering ~30,000 cards for genes with no signal anywhere took tens of minutes
and buried the reader; the first real OsHV-1 pool renders a few thousand
instead. Pass `--all-cards` to render every candidate regardless.

Each card shows the tier, the contributing studies, the pooled and
BH-adjusted p-values with the family size, and a forest-style table of
per-study effect estimates with each study's own adjusted p, so a reviewer can
see exactly which inputs drove the score rather than trusting a single number.

## Related documentation

- [interpreting_meta_analysis.md](interpreting_meta_analysis.md)
- [resilience_vs_exposure.md](resilience_vs_exposure.md)
- [identifier_mapping.md](identifier_mapping.md)
- [design.md](design.md#6-candidate-scoring-framework-transparent-not-black-box)
