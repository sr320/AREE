# Resilience vs. exposure-only evidence

A recurring failure mode in biomarker discovery is treating "this gene's
expression changed when we stressed the animal" as if it were equivalent to
"this gene predicts which animals will survive the stress." AREE refuses to
make that substitution silently. This page explains the mechanism.

## The four `resilience_classification` categories

Every comparison in a study record (`comparisons[].resilience_classification`
in `schemas/study.schema.json`) is tagged with exactly one of:

- **`resilience`** — the phenotype is a direct measure of organismal capacity
  to survive or perform under adversity (survival, thermal tolerance, growth
  under stress, larval viability, salinity/hypoxia/acidification tolerance,
  recovery following stress, metabolic resilience, reproductive performance,
  disease resistance).
- **`stress_response`** — the measurement documents a molecular or
  physiological response to a stressor, but is not itself an outcome measure
  of resilience (e.g. immune responsiveness on its own, without a linked
  survival/performance outcome).
- **`disease`** — a disease susceptibility/resistance/pathogen-load outcome.
  Related to resilience but kept as its own category because disease
  challenge studies have their own confounds (dose, pathogen strain,
  challenge route) distinct from abiotic stress.
- **`exposure_only`** — the record describes an exposure or treatment
  condition, not an outcome phenotype at all.

This same four-way split is baked into the phenotype ontology itself via each
term's `resilience_relevance` field (`registry/controlled_vocabularies/phenotype_ontology.yaml`)
— see [defining_a_phenotype.md](defining_a_phenotype.md). A comparison's
`resilience_classification` should normally match its phenotype term's
`resilience_relevance`; divergence should be justified in curation notes.

## Why this distinction is enforced, not advisory

Design document [design.md](design.md), sections 2 and 6, treats this as a
structural requirement rather than a style guideline:

- Section 2 (data model) locates phenotype and stressor fields at the
  **comparison** level specifically so a study cannot be summarized by a
  single ambiguous label — "changed under stress" and "predicted survival"
  never get flattened into the same field.
- Section 6 (candidate scoring) uses `resilience_relevance` directly as an
  input to `phenotype_relevance_score`
  (`src/prioritize/scoring.py`, `PHENOTYPE_RELEVANCE_SCORE`):

  | relevance | score |
  |---|---|
  | `resilience` | 1.0 |
  | `disease` | 0.8 |
  | `stress_response` | 0.4 |
  | `exposure_only` | 0.1 |

  A gene that only ever shows up in `stress_response` or `exposure_only`
  comparisons is structurally capped at a low `phenotype_relevance_score` no
  matter how large its effect size or how many studies replicate it. It
  cannot reach the `high_priority_cross_study` tier purely on statistical
  strength — the tiering gate in `src/prioritize/rank.py` requires
  `phenotype_relevance_score > 0.1`, which `exposure_only` evidence alone
  does not clear. See
  [interpreting_candidate_scores.md](interpreting_candidate_scores.md) for
  the full gating logic.

## What this means in practice when curating a study

- Do not default every "differentially expressed under stressor X" result to
  `resilience_classification: resilience` just because the study is about a
  stressor. Ask: does this comparison's phenotype measure organismal
  capacity to survive/perform (resilience), a downstream molecular response
  (stress_response), or a disease outcome (disease)? If the "phenotype" is
  actually just describing the treatment condition itself, it is
  `exposure_only`.
- A single study can legitimately contain comparisons in more than one
  category — e.g. `GIGAS_PATH03` (`registry/studies/GIGAS_PATH03.yaml`) is
  classified `disease` because its phenotype is `disease_susceptibility`
  measured via tissue pathology score rather than a survival curve, which is
  explicitly noted in its `limitations` field.
- The point is never to discard `stress_response` or `exposure_only`
  evidence — it stays in the evidence table and is visible in evidence cards
  — but AREE will not let it masquerade as validated resilience evidence.

## Related documentation

- [defining_a_phenotype.md](defining_a_phenotype.md)
- [interpreting_candidate_scores.md](interpreting_candidate_scores.md)
- [design.md](design.md)
