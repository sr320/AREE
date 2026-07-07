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
    "assay_diversity_score": 10,     # multi-omics convergence
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


def compute_components(meta_row: dict, mapping_confidences: list, quality_flags: list) -> dict:
    """Compute the 0-1 component scores for one candidate (one meta-analysis row).

    meta_row is expected to have the columns produced by
    meta_analysis.run.run_meta_analysis (k_studies, total_sample_size,
    pooled_effect, p_value, direction_consistency, distinct_tissues,
    distinct_life_stages, i_squared, phenotype).
    """
    phenotype_vocab = load_vocab("phenotype_ontology")
    relevance = phenotype_vocab.get(meta_row["phenotype"], {}).get("resilience_relevance", "exposure_only")

    p = max(meta_row["p_value"], 1e-300)
    significance_score = _clip01(-math.log10(p) / 5.0)  # p=1e-5 -> 1.0

    worst_mapping = min(
        (MAPPING_CONFIDENCE_SCORE.get(m, 0.0) for m in mapping_confidences), default=0.0
    )

    return {
        "n_studies_score": _clip01(meta_row["k_studies"] / 5.0),
        "sample_size_score": _clip01(meta_row["total_sample_size"] / 100.0),
        "effect_magnitude_score": _clip01(abs(meta_row["pooled_effect"]) / 2.0),
        "significance_score": significance_score,
        "direction_consistency_score": _clip01(meta_row["direction_consistency"]),
        "phenotype_relevance_score": PHENOTYPE_RELEVANCE_SCORE.get(relevance, 0.1),
        "context_breadth_score": _clip01(
            (meta_row["distinct_tissues"] * meta_row["distinct_life_stages"]) / 3.0
        ),
        "assay_diversity_score": _clip01(meta_row.get("n_distinct_assays", 1) / 3.0),
        "mapping_confidence_score": worst_mapping,
        "quality_score": _clip01(1 - len(quality_flags) / 5.0),
        "heterogeneity_penalty": _clip01(meta_row.get("i_squared", 0.0) / 100.0),
    }


def candidate_score(components: dict) -> float:
    """Pure function: same components in, same score out. 0-100 scale."""
    positive = sum(WEIGHTS[name] * components[name] for name in WEIGHTS)
    penalty = HETEROGENEITY_PENALTY_WEIGHT * components.get("heterogeneity_penalty", 0.0)
    return round(max(0.0, min(100.0, positive - penalty)), 2)
