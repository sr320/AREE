"""Cross-species candidates: orthogroups whose evidence converges across species. WP6.

Input is the orthogroup-level meta-analysis (meta_analysis/orthogroups.py): one
row per orthogroup, phenotype and feature type, pooling each species' own
gene-level estimate. This module tiers and scores those rows the same way
prioritize/rank.py does for single-species candidates: **hard gates first,
score second**. A high score never promotes a row past a gate it fails.

Tiers, strongest first:

``cross_species_convergence``
    Every gate passes: the WP6 replication gate (at least two species, each
    with k >= 2 studies of its own), every species' estimate has the same sign,
    the pooled cross-species estimate survives BH adjustment (q <= 0.05), the
    weakest contributing evidence class is more than exposure-only, and the
    quality flags are within the gene-level limit.
``cross_species_emerging``
    Species agree in sign and the pooled estimate is significant, but at least
    one other gate fails. Usually the replication gate: one species rests on a
    single study. ``gates_failed`` names each failed gate.
``cross_species_discordant``
    The species' estimates disagree in sign. Surfaced, not dropped: an ortholog
    responding in opposite directions in two species is a finding about the
    biology, or about the studies, and either way a reader must see it.
``cross_species_unsupported``
    Species agree in sign but the pooled estimate is not significant.

Agreement in sign is judged on point estimates, as at gene level, so a
near-zero estimate of opposite sign makes a row discordant. The card shows
each species' estimate and confidence interval so a reader can tell a real
reversal from noise around zero.
"""
from __future__ import annotations

import math

import pandas as pd

from .rank import EVIDENCE_CLASS_BY_RELEVANCE, HIGH_PRIORITY_MIN_QUALITY_SCORE, _partition_value
from .scoring import (
    MAPPING_CONFIDENCE_SCORE,
    PHENOTYPE_RELEVANCE_SCORE,
    _clip01,
    _saturating,
    conservative_relevance,
    significance_p_value,
)

TIER_ORDER = [
    "cross_species_convergence",
    "cross_species_emerging",
    "cross_species_discordant",
    "cross_species_unsupported",
]
CONVERGENCE_MAX_ADJUSTED_P = 0.05
CONVERGENCE_MIN_RELEVANCE_SCORE = PHENOTYPE_RELEVANCE_SCORE["exposure_only"]  # must exceed

# Weights sum to 100; between-species heterogeneity is a separate penalty.
WEIGHTS = {
    "n_species_score": 15,                # more species agreeing is stronger conservation evidence
    "per_species_replication_score": 20,  # the weakest species' study count: a chain is its weakest link
    "sample_size_score": 5,               # total biological replication across species
    "effect_magnitude_score": 10,         # pooled cross-species effect
    "significance_score": 10,             # BH-adjusted pooled p, within its family
    "direction_consistency_score": 15,    # share of species agreeing with the majority sign
    "phenotype_relevance_score": 10,      # weakest evidence class among all contributing records
    "orthology_score": 8,                 # how like-for-like the orthogroup is
    "mapping_confidence_score": 4,        # worst identifier mapping behind any species' estimate
    "quality_score": 3,                   # quality flags across all contributing records
}
HETEROGENEITY_PENALTY_WEIGHT = 15  # subtracted, scaled by between-species I^2 / 100

# Single-copy in every species is the cleanest one-to-one comparison; single-copy
# among the species present is nearly as good; the demo crosswalk is one gene
# per species by construction; multi-copy (opt-in only) compares one paralog
# with another species' different paralog.
ORTHOLOGY_SCORE = {
    "single_copy_all": 1.0,
    "single_copy": 0.9,
    "demo_crosswalk": 0.9,
    "multi_copy": 0.3,
}

GATES = ["replication", "direction", "significance", "phenotype_relevance", "quality"]


def _split(value) -> list[str]:
    if value is None or (isinstance(value, float) and value != value):
        return []
    return [x for x in str(value).split("|") if x]


def species_estimates_for(species_df: pd.DataFrame, row: dict) -> pd.DataFrame:
    """The per-species estimates behind one orthogroup result row."""
    sim = _partition_value(row["simulated"])
    mask = (
        (species_df["orthogroup_id"].astype(str) == str(row["orthogroup_id"]))
        & (species_df["phenotype"] == row["phenotype"])
        & (species_df["feature_type"] == row["feature_type"])
        & (species_df["simulated"].map(_partition_value) == sim)
    )
    return species_df[mask]


def compute_components(row: dict, species: pd.DataFrame, relevance: str) -> dict:
    """0-1 component scores for one orthogroup row. Pure: same inputs, same output."""
    k = pd.to_numeric(species["k_studies"], errors="coerce").fillna(0)
    mapping = {m for v in species.get("mapping_confidences", []) for m in _split(v)}
    flags = {f for v in species.get("quality_flags_union", []) for f in _split(v)}
    n = pd.to_numeric(species.get("total_sample_size", pd.Series(dtype=float)), errors="coerce").fillna(0)
    i2 = row.get("i_squared_between_species")
    return {
        "n_species_score": _saturating(row["n_species"], 2.0),
        "per_species_replication_score": _saturating(k.min() if len(k) else 0, 2.0),
        "sample_size_score": _saturating(n.sum(), 100.0),
        "effect_magnitude_score": _saturating(abs(row["pooled_effect"]), 2.0),
        "significance_score": _saturating(-math.log10(significance_p_value(row)), 5.0),
        "direction_consistency_score": _clip01(row["direction_consistency_species"]),
        "phenotype_relevance_score": PHENOTYPE_RELEVANCE_SCORE.get(relevance, 0.1),
        "orthology_score": ORTHOLOGY_SCORE.get(row["copy_class"], 0.0),
        "mapping_confidence_score": min((MAPPING_CONFIDENCE_SCORE.get(m, 0.0) for m in mapping), default=0.0),
        "quality_score": _clip01(1 - len(flags) / 5.0),
        "heterogeneity_penalty": _clip01((i2 if i2 == i2 and i2 is not None else 0.0) / 100.0),
    }


def cross_species_score(components: dict) -> float:
    positive = sum(WEIGHTS[name] * components[name] for name in WEIGHTS)
    penalty = HETEROGENEITY_PENALTY_WEIGHT * components.get("heterogeneity_penalty", 0.0)
    return round(max(0.0, min(100.0, positive - penalty)), 2)


def gates_failed(row: dict, components: dict) -> list[str]:
    """Each convergence gate this row fails, in GATES order."""
    failed = []
    if not bool(row["meets_replication_gate"]):
        failed.append("replication")
    if row["direction_consistency_species"] < 1.0:
        failed.append("direction")
    if significance_p_value(row) > CONVERGENCE_MAX_ADJUSTED_P:
        failed.append("significance")
    if components["phenotype_relevance_score"] <= CONVERGENCE_MIN_RELEVANCE_SCORE:
        failed.append("phenotype_relevance")
    if components["quality_score"] < HIGH_PRIORITY_MIN_QUALITY_SCORE:
        failed.append("quality")
    return failed


def tier_for(failed: list[str]) -> str:
    if not failed:
        return "cross_species_convergence"
    if "direction" in failed:
        return "cross_species_discordant"
    if "significance" in failed:
        return "cross_species_unsupported"
    return "cross_species_emerging"


def _relevance(row: dict, species: pd.DataFrame, classification_by_id: dict) -> str:
    ids = [i for v in species.get("contributing_evidence_ids", []) for i in _split(v)]
    return conservative_relevance(row["phenotype"], [classification_by_id.get(i) for i in ids])


def build_cross_species_candidates(
    og_results: pd.DataFrame, species_df: pd.DataFrame, evidence: pd.DataFrame
) -> pd.DataFrame:
    """Tier and score every orthogroup result row. Every row is returned."""
    if og_results is None or len(og_results) == 0:
        return pd.DataFrame()
    classification_by_id = {}
    if "resilience_classification" in evidence.columns:
        classification_by_id = dict(zip(evidence["evidence_id"].astype(str), evidence["resilience_classification"]))

    rows = []
    for row in og_results.to_dict("records"):
        species = species_estimates_for(species_df, row)
        relevance = _relevance(row, species, classification_by_id)
        components = compute_components(row, species, relevance)
        failed = gates_failed(row, components)
        rows.append({
            **row,
            "tier": tier_for(failed),
            "score": cross_species_score(components),
            "gates_failed": "|".join(failed),
            "evidence_class": EVIDENCE_CLASS_BY_RELEVANCE.get(relevance, "exposure_only"),
            "min_k_studies_per_species": int(pd.to_numeric(species["k_studies"]).min()) if len(species) else 0,
            "total_sample_size": int(pd.to_numeric(species["total_sample_size"], errors="coerce").fillna(0).sum())
            if "total_sample_size" in species.columns else 0,
            "mapping_confidences": "|".join(sorted({m for v in species.get("mapping_confidences", [])
                                                    for m in _split(v)})),
            "quality_flags_union": "|".join(sorted({f for v in species.get("quality_flags_union", [])
                                                    for f in _split(v)})),
            **{f"component_{k}": v for k, v in components.items()},
        })
    out = pd.DataFrame(rows)
    rank = out["tier"].map({t: i for i, t in enumerate(TIER_ORDER)})
    return (out.assign(_rank=rank).sort_values(["_rank", "score"], ascending=[True, False])
            .drop(columns="_rank").reset_index(drop=True))
