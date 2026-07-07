# Data governance and provenance

AREE's evidentiary value depends on every number being traceable back to a
source file and a set of parameters. This page describes what provenance
information is mandatory, how missing information is handled, how curation
decisions are logged, and the causality boundary the resource enforces.

## Mandatory provenance fields

Per [design.md](design.md#4-provenance-model), every workflow output and
every evidence record carries:

- source accession(s) and `study_id`
- input file name(s) and checksum(s) — evidence records store this as
  `source_file`, a relative path plus a sha256 checksum
  (`schemas/evidence.schema.json` documents the convention
  `'path#sha256:...'`; `src/common/sha256_file` computes it)
- parameter set (explicit key/value, not "defaults")
- workflow name and semantic version (`workflow_version` on every evidence
  record)
- tool versions, from the container manifest or environment lock
- reference genome/annotation version (`genome_assembly`,
  `annotation_version`)
- `date_generated` (ISO 8601) — supplied explicitly by the caller. The
  codebase does not call wall-clock time internally for this field, so
  provenance timestamps are always explicit and reproducible rather than
  drifting with when a command happens to be re-run
- `generated_by` — either `automated:<workflow_id>@<version>` or
  `curator:<name>`
- any manual curation decisions, recorded as a separate append-only log
  entry rather than an in-place edit

At the study level, `provenance.registered_by` and `provenance.date_registered`
are required by `schemas/study.schema.json`; `provenance.curation_notes` and
`provenance.source_links` are available for anything that needs a paper
trail beyond the structured fields.

## Missingness is explicit, never silent

If a source study does not report a field, that field is recorded as `null`
(with a `missing_reason` where feasible), not omitted from the schema and not
imputed. Several schema fields are typed `["string", "null"]` or similar
specifically to make "not reported" a representable, queryable state rather
than an absence that looks the same as "not yet entered." Quality flags
(`registry/controlled_vocabularies/quality_flags.yaml`) exist for the same
reason at a coarser grain — `processed_only`, `low_replication`,
`unbalanced_design`, `no_multiple_testing_correction`,
`ambiguous_phenotype_definition`, `single_timepoint`,
`identifier_mapping_uncertain` are all first-class, visible labels rather
than caveats buried in free text.

Treatment and stressor information follows the same principle at a different
angle: the exact original wording (`stressor_original`) is always preserved
alongside the standardized ontology mapping (`stressor_standardized`) — see
[identifier_mapping.md](identifier_mapping.md) and
[defining_a_phenotype.md](defining_a_phenotype.md) for the identifier and
phenotype analogs of this same original-vs-standardized pattern. Mapping is
additive, never destructive; nothing overwrites the source study's own
description.

## Curation decisions are logged, not silently applied

Where a curator makes a judgment call — for example, classifying a
comparison's `resilience_classification` differently from its phenotype
term's default `resilience_relevance` (see
[resilience_vs_exposure.md](resilience_vs_exposure.md)), or resolving an
ambiguous identifier via `data/mappings/ambiguous_symbol_map.yaml` — that
decision is recorded in `provenance.curation_notes` or `limitations` rather
than applied invisibly. The goal is that a second curator or reviewer can
reconstruct why a record looks the way it does without re-deriving the
judgment from scratch.

## Contradictory findings are surfaced, not removed

AREE does not delete or suppress a study's results because they conflict
with another study's. The conflicting-evidence example documented in
[interpreting_meta_analysis.md](interpreting_meta_analysis.md) (`sod1` /
`LOC105331241` showing opposite-direction effects in `GIGAS_OA02` and
`GIGAS_LARV05`) is retained in full in the evidence table, is visible in the
pooled meta-analysis output via `direction_consistency` and `i_squared`, and
would appear in an evidence card if that candidate reached scoring — labeled
as low-consistency, not hidden or averaged away without explanation.

## The causality boundary

AREE identifies **associations and evidence convergence**, not confirmed
mechanistic causation. A candidate with strong, replicated, direction-consistent
support across studies is a well-supported candidate for further validation —
not a validated biomarker. This is enforced structurally, not just stated as
a caveat:

- No candidate reaches the `high_priority_cross_study` tier without meeting
  explicit replication and consistency gates (see
  [interpreting_candidate_scores.md](interpreting_candidate_scores.md)), and
  even that tier's evidence-card language and recommended next step
  (`src/reporting/evidence_cards.py`, `NEXT_STEP_BY_TIER`) point toward
  further validation, not confirmation.
- `emerging` candidates — single-study support, however strong — are always
  explicitly labeled as requiring replication and are structurally excluded
  from the higher tiers regardless of score.
- Statistical significance in one study is never, by construction, sufficient
  on its own to reach the top tier.

## Related documentation

- [interpreting_candidate_scores.md](interpreting_candidate_scores.md)
- [interpreting_meta_analysis.md](interpreting_meta_analysis.md)
- [resilience_vs_exposure.md](resilience_vs_exposure.md)
- [design.md](design.md#4-provenance-model)
