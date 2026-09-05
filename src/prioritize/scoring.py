"""Transparent candidate scoring.

Every component below is a named, independently inspectable 0-1 value computed
from fields already present in the meta-analysis/evidence tables. The final
score is a plain weighted sum minus a heterogeneity penalty — there is no
hidden model. Adjust WEIGHTS (and leave a comment explaining why) rather than
special-casing individual candidates.

This module deliberately does NOT decide tier membership (high-priority vs
multi-omics vs emerging) — see prioritize/rank.py for the hard gating rules
that a high score alone cannot override.
"""
from __future__ import annotations

import math
from functools import lru_cache

from common import load_vocab

# Weights sum to 100 across the positive components; heterogeneity is a
# separate penalty subtracted from the total, not part of the 100.
WEIGHTS = {
    "n_studies_score": 20,           # independent replication is the single strongest signal
    "sample_size_score": 10,         # total biological replication behind the pooled estimate
    "effect_magnitude_score": 10,    # larger effects are more actionable, but not sufficient alone
    "significance_score": 10,        # adjusted-significance strength
    "direction_consistency_score": 20,  # agreement in direction across studies matters as much as study count
    "phenotype_relevance_score": 10, # resilience > disease > stress_response > exposure_only
    "context_breadth_score": 5,      # spread across tissues/life stages
    "assay_diversity_score": 10,     # molecular layers with significant same-phenotype support
    "mapping_confidence_score": 3,   # trust in the identifier harmonization itself
    "quality_score": 2,              # study/data quality flags
}
HETEROGENEITY_PENALTY_WEIGHT = 15  # subtracted, scaled by I^2 / 100

PHENOTYPE_RELEVANCE_SCORE = {
    "resilience": 1.0,
    "disease": 0.8,
    "stress_response": 0.4,
    "exposure_only": 0.1,
}

MAPPING_CONFIDENCE_SCORE = {
    "exact": 1.0,
    "one_to_one_ortholog": 0.9,
    "one_to_many_ortholog": 0.5,
    "many_to_one_ortholog": 0.5,
    "inferred": 0.6,
    "unresolved": 0.0,
}


def _clip01(x: float) -> float:
    return max(0.0, min(1.0, x))


def _saturating(value: float, half: float) -> float:
    """Map an unbounded non-negative quantity into [0, 1) without a ceiling.

    Returns `v / (v + half)`: strictly increasing over the whole domain, equal
    to 0.5 at `half`, and asymptotic to — never equal to — 1.0.

    This replaces `min(1, v / ceiling)` for every component whose input is
    unbounded. Clipping was fine on the simulated demo, where nothing reached
    the ceiling. On the first genome-wide real pool it broke the top of the
    ranking: `|pooled_effect|` clipped at 2 and adjusted p clipped at 1e-5, so
    126 of the 1,994 high-priority candidates scored an identical 68.57 despite
    effects spanning 2 to 9.7 log2FC and adjusted p spanning 1e-5 to 1e-75 —
    ties in exactly the region a reader is choosing validation targets from.

    Each `half` below is the value that used to be the ceiling, so the ordering
    of the old score is preserved and only the saturation is removed. Scores are
    therefore lower in absolute terms than before 2026-09-05; compare candidates
    within a run, not across this change.
    """
    v = max(0.0, float(value))
    return v / (v + half)


@lru_cache(maxsize=1)
def _phenotype_relevance() -> dict:
    """phenotype term id -> resilience_relevance class, loaded once.

    Scoring is called once per candidate, and a genome-wide pool has tens of
    thousands of candidates; re-reading the ontology YAML on each call was a
    measurable share of ranking time.
    """
    return {
        term_id: term.get("resilience_relevance", "exposure_only")
        for term_id, term in load_vocab("phenotype_ontology").items()
    }


def phenotype_relevance(phenotype: str) -> str:
    """The phenotype's `resilience_relevance` class, defaulting to exposure_only."""
    return _phenotype_relevance().get(phenotype, "exposure_only")


def significance_p_value(meta_row: dict) -> float:
    """The p-value the score and the tier gates judge significance on.

    The BH-adjusted pooled p is preferred; the raw pooled p is the fallback for
    callers that build a meta row by hand without running the family-wise
    adjustment. Either way the result is bounded away from zero so log10 is
    defined.
    """
    adjusted = meta_row.get("adjusted_p_value")
    p = adjusted if adjusted is not None and adjusted == adjusted else meta_row["p_value"]
    return max(float(p), 1e-300)


def compute_components(meta_row: dict, mapping_confidences: list, quality_flags: list) -> dict:
    """Compute the 0-1 component scores for one candidate (one meta-analysis row).

    meta_row is expected to have the columns produced by
    meta_analysis.run.run_meta_analysis (k_studies, total_sample_size,
    pooled_effect, p_value, adjusted_p_value, direction_consistency,
    distinct_tissues, distinct_life_stages, i_squared, phenotype) plus
    n_supporting_layers, which prioritize.rank adds.
    """
    relevance = phenotype_relevance(meta_row["phenotype"])

    # -log10(q) is itself unbounded, so saturate on that scale: q=1e-5 -> 0.5.
    significance_score = _saturating(-math.log10(significance_p_value(meta_row)), 5.0)

    worst_mapping = min(
        (MAPPING_CONFIDENCE_SCORE.get(m, 0.0) for m in mapping_confidences), default=0.0
    )

    return {
        # Unbounded inputs: saturating, so more studies / bigger n / bigger
        # effect / smaller q always ranks strictly higher.
        "n_studies_score": _saturating(meta_row["k_studies"], 5.0),
        "sample_size_score": _saturating(meta_row["total_sample_size"], 100.0),
        "effect_magnitude_score": _saturating(abs(meta_row["pooled_effect"]), 2.0),
        "significance_score": significance_score,
        "context_breadth_score": _saturating(
            meta_row["distinct_tissues"] * meta_row["distinct_life_stages"], 3.0
        ),
        # Bounded inputs keep the clip: the number of molecular layers is capped
        # by the feature_types vocabulary, and the rest are already 0-1 by
        # construction, so clipping them saturates nothing that could go higher.
        "assay_diversity_score": _clip01(meta_row.get("n_supporting_layers", 1) / 3.0),
        "direction_consistency_score": _clip01(meta_row["direction_consistency"]),
        "phenotype_relevance_score": PHENOTYPE_RELEVANCE_SCORE.get(relevance, 0.1),
        "mapping_confidence_score": worst_mapping,
        "quality_score": _clip01(1 - len(quality_flags) / 5.0),
        "heterogeneity_penalty": _clip01(meta_row.get("i_squared", 0.0) / 100.0),
    }


def candidate_score(components: dict) -> float:
    """Pure function: same components in, same score out. 0-100 scale."""
    positive = sum(WEIGHTS[name] * components[name] for name in WEIGHTS)
    penalty = HETEROGENEITY_PENALTY_WEIGHT * components.get("heterogeneity_penalty", 0.0)
    return round(max(0.0, min(100.0, positive - penalty)), 2)
