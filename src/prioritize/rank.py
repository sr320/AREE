"""Candidate tiering: hard gates first, score second.

A candidate's numeric score never promotes it past these gates — this is the
mechanism that keeps "significant in one study" from reading as "validated"
(see docs/design.md section 6).
"""
from __future__ import annotations

import ast

import pandas as pd

from .scoring import candidate_score, compute_components

HIGH_PRIORITY_MIN_STUDIES = 2
HIGH_PRIORITY_MIN_DIRECTION_CONSISTENCY = 0.7
HIGH_PRIORITY_MIN_QUALITY_SCORE = 0.4
MULTI_OMICS_MIN_ASSAYS = 2


def _parse_flags(value) -> list:
    if isinstance(value, list):
        return value
    if isinstance(value, str) and value.startswith("["):
        return ast.literal_eval(value)
    return []


def _assay_diversity_map(evidence_df: pd.DataFrame) -> dict:
    mapped = evidence_df[evidence_df["feature_id_standardized"].notna()]
    return mapped.groupby("feature_id_standardized")["feature_type"].nunique().to_dict()


def _evidence_subset(evidence_df: pd.DataFrame, feature_id: str, phenotype: str, feature_type: str) -> pd.DataFrame:
    return evidence_df[
        (evidence_df["feature_id_standardized"] == feature_id)
        & (evidence_df["phenotype"] == phenotype)
        & (evidence_df["feature_type"] == feature_type)
    ]


def build_candidates(meta_df: pd.DataFrame, evidence_df: pd.DataFrame) -> pd.DataFrame:
    if len(meta_df) == 0:
        return pd.DataFrame()

    assay_diversity = _assay_diversity_map(evidence_df)

    rows = []
    for _, meta_row in meta_df.iterrows():
        meta_dict = meta_row.to_dict()
        subset = _evidence_subset(
            evidence_df, meta_dict["feature_id_standardized"], meta_dict["phenotype"], meta_dict["feature_type"]
        )
        mapping_confidences = sorted(subset["mapping_confidence"].unique())
        quality_flags = sorted({f for flags in subset["quality_flags"] for f in _parse_flags(flags)})
        n_distinct_assays = assay_diversity.get(meta_dict["feature_id_standardized"], 1)

        meta_dict["n_distinct_assays"] = n_distinct_assays
        components = compute_components(meta_dict, mapping_confidences, quality_flags)
        score = candidate_score(components)

        is_high_priority = (
            meta_dict["k_studies"] >= HIGH_PRIORITY_MIN_STUDIES
            and components["phenotype_relevance_score"] > 0.1
            and meta_dict["direction_consistency"] >= HIGH_PRIORITY_MIN_DIRECTION_CONSISTENCY
            and components["quality_score"] >= HIGH_PRIORITY_MIN_QUALITY_SCORE
        )
        is_multi_omics = n_distinct_assays >= MULTI_OMICS_MIN_ASSAYS

        if is_high_priority:
            tier = "high_priority_cross_study"
        elif is_multi_omics:
            tier = "multi_omics_convergence"
        else:
            tier = "emerging"

        rows.append({
            **meta_dict,
            "score": score,
            "tier": tier,
            "is_high_priority": is_high_priority,
            "is_multi_omics_convergence": is_multi_omics,
            "n_distinct_assays": n_distinct_assays,
            "mapping_confidences": "|".join(mapping_confidences),
            "quality_flags_union": "|".join(quality_flags),
            **{f"component_{k}": v for k, v in components.items()},
        })

    out = pd.DataFrame(rows)
    return out.sort_values(["tier", "score"], ascending=[True, False])
