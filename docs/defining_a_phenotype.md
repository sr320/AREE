# Defining a phenotype

AREE uses a controlled phenotype ontology so that "survival," "mortality,"
and "thermal tolerance" mean the same thing across every study that reports
them, instead of each study inventing its own free-text label.

## Where it lives

`registry/controlled_vocabularies/phenotype_ontology.yaml`. Each term has:

```yaml
- id: thermal_tolerance
  label: Thermal tolerance
  definition: Capacity to maintain function/survival under elevated or fluctuating temperature.
  resilience_relevance: resilience
  unit_examples: ["CTmax_C", "percent_survival_at_temp"]
```

- `id` — the stable identifier used in `comparisons[].phenotype` and in the
  evidence table's `phenotype` column. This is what code and schemas
  reference, not the human-readable `label`.
- `label` — display name.
- `definition` — a one-sentence operational definition, so curators applying
  the term to a new study can judge whether it fits.
- `resilience_relevance` — one of `resilience`, `stress_response`, `disease`,
  `exposure_only`. This is the field that keeps AREE from treating "this gene
  changed under stress" as "this gene predicts resilience." See
  [resilience_vs_exposure.md](resilience_vs_exposure.md) for the full
  rationale.
- `unit_examples` — illustrative, not enforced, units a study might report
  this phenotype in (`phenotype_unit` on the comparison is free text).

The current ontology (version `0.1.0`) has 15 terms: `survival`, `mortality`,
`thermal_tolerance`, `growth_under_stress`, `pathogen_load`,
`disease_susceptibility`, `disease_resistance`, `reproductive_performance`,
`larval_viability`, `metabolic_resilience`, `recovery_following_stress`,
`immune_responsiveness`, `acidification_tolerance`, `salinity_tolerance`,
`hypoxia_tolerance`.

## Using an existing term

In a study's `comparisons[]` entry, set `phenotype` to the term `id` (not the
label), e.g. `phenotype: "larval_viability"`. The value is checked against
this vocabulary by `aree validate-study`.

`resilience_classification` on the same comparison is a separate field (see
`schemas/study.schema.json`) and should normally agree with the ontology
term's `resilience_relevance` — the schema comment notes it "must agree...
unless explicitly overridden with a note." If you deliberately set a
comparison's classification differently from the ontology default (for
example, treating a normally-`resilience` phenotype as `exposure_only` for a
particular under-powered comparison), record why in `limitations` or
`provenance.curation_notes`.

## Proposing a new term

If a study reports a resilience-relevant outcome that doesn't fit any
existing term:

1. Add a new entry to `registry/controlled_vocabularies/phenotype_ontology.yaml`
   with a unique `id` (snake_case), a `label`, a precise `definition`, and a
   `resilience_relevance` classification. Be conservative about
   `resilience_relevance` — see
   [resilience_vs_exposure.md](resilience_vs_exposure.md) before defaulting
   to `resilience`.
2. Bump the vocabulary file's `version` field.
3. Reference the new term `id` in your study's `comparisons[].phenotype`.
4. Run `aree validate-study` on any studies using the new term to confirm it
   resolves.
5. Note the addition in your PR description so reviewers can check the
   `resilience_relevance` call — this is a scientific judgment, not a
   mechanical one, and it directly affects candidate scoring (see
   [interpreting_candidate_scores.md](interpreting_candidate_scores.md),
   `phenotype_relevance_score`).

Avoid creating near-duplicate terms (e.g. a second "heat tolerance" term)
when an existing term's `definition` already covers the case — prefer reusing
`thermal_tolerance` and letting `stressor_original`/`exposure_*` fields carry
the specifics.

## Related documentation

- [resilience_vs_exposure.md](resilience_vs_exposure.md)
- [adding_a_study.md](adding_a_study.md)
- [design.md](design.md#3-controlled-vocabularies)
