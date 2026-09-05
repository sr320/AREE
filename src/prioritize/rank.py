"""Candidate tiering: hard gates first, score second.

A candidate's numeric score never promotes it past these gates — this is the
mechanism that keeps "significant in one study" from reading as "validated"
(see docs/design.md section 6).

Ranking is a single grouped pass over the evidence table. Each candidate's
evidence subset is looked up from a precomputed index instead of re-filtering
the whole table per candidate, which is what makes a genome-wide pool of
30,000+ candidates rank in seconds rather than tens of minutes.
"""
from __future__ import annotations

import ast

import pandas as pd

from .scoring import candidate_score, compute_components, significance_p_value

HIGH_PRIORITY_MIN_STUDIES = 2
HIGH_PRIORITY_MIN_DIRECTION_CONSISTENCY = 0.7
HIGH_PRIORITY_MIN_QUALITY_SCORE = 0.4
# Two genome-wide studies agree in sign for roughly half of all null genes, so
# study count and direction consistency alone cannot define the top tier. The
# pooled effect must also survive family-wise false-discovery control.
HIGH_PRIORITY_MAX_ADJUSTED_P = 0.05
MULTI_OMICS_MIN_ASSAYS = 2

# Order in which tiers are presented: strongest evidence first.
TIER_ORDER = ["high_priority_cross_study", "multi_omics_convergence", "emerging"]

PARTITION_COLUMNS = ["_simulated", "_species_taxid"]
CANDIDATE_GROUP_COLUMNS = ["feature_id_standardized", "phenotype", "feature_type", *PARTITION_COLUMNS]
FEATURE_GROUP_COLUMNS = ["feature_id_standardized", *PARTITION_COLUMNS]


def _parse_flags(value) -> list:
    if isinstance(value, list):
        return value
    if isinstance(value, str) and value.startswith("["):
        return ast.literal_eval(value)
    return []


def _partition_value(value) -> str:
    if isinstance(value, bool):
        return str(value).lower()
    return str(value).strip().lower()


def add_partition_columns(evidence_df: pd.DataFrame) -> pd.DataFrame:
    """Return the evidence with normalized `_simulated` / `_species_taxid` columns.

    `simulated` arrives as a bool in memory but as the text "True"/"False" from a
    TSV, and `species_taxid` as int or text; normalizing once here is what lets
    every downstream grouping key be compared by plain equality. Idempotent.
    """
    if all(col in evidence_df.columns for col in PARTITION_COLUMNS):
        return evidence_df
    return evidence_df.assign(
        _simulated=evidence_df["simulated"].map(_partition_value),
        _species_taxid=evidence_df["species_taxid"].map(_partition_value),
    )


def candidate_key(meta_row: dict) -> tuple:
    return (
        str(meta_row["feature_id_standardized"]),
        meta_row["phenotype"],
        meta_row["feature_type"],
        _partition_value(meta_row["simulated"]),
        _partition_value(meta_row["species_taxid"]),
    )


class EvidenceIndex:
    """Positional indices into a normalized evidence frame, built once per run."""

    def __init__(self, evidence_df: pd.DataFrame):
        evidence = add_partition_columns(evidence_df)
        mapped = evidence[evidence["feature_id_standardized"].notna()].copy()
        mapped["feature_id_standardized"] = mapped["feature_id_standardized"].astype(str)
        self.evidence = mapped
        self._by_candidate = mapped.groupby(CANDIDATE_GROUP_COLUMNS, sort=False).indices
        self._by_feature = mapped.groupby(FEATURE_GROUP_COLUMNS, sort=False).indices
        self.assay_diversity = (
            mapped.groupby(FEATURE_GROUP_COLUMNS)["feature_type"].nunique().to_dict()
        )

    def _slice(self, index, key) -> pd.DataFrame:
        positions = index.get(key)
        if positions is None:
            return self.evidence.iloc[0:0]
        return self.evidence.iloc[positions]

    def for_candidate(self, meta_row: dict) -> pd.DataFrame:
        """All records for this feature / phenotype / feature type / partition."""
        return self._slice(self._by_candidate, candidate_key(meta_row))

    def contributing(self, meta_row: dict) -> pd.DataFrame:
        """The records that actually entered the pooled estimate."""
        subset = self.for_candidate(meta_row)
        ids = meta_row.get("contributing_evidence_ids")
        if not ids or ids != ids:
            return subset
        wanted = set(str(ids).split("|"))
        return subset[subset["evidence_id"].astype(str).isin(wanted)]

    def other_feature_types(self, meta_row: dict) -> pd.DataFrame:
        """Records for the same feature and partition from OTHER feature types."""
        key = candidate_key(meta_row)
        same_feature = self._slice(self._by_feature, (key[0], key[3], key[4]))
        return same_feature[same_feature["feature_type"] != meta_row["feature_type"]]

    def n_distinct_assays(self, meta_row: dict) -> int:
        key = candidate_key(meta_row)
        return int(self.assay_diversity.get((key[0], key[3], key[4]), 1))


def build_candidates(meta_df: pd.DataFrame, evidence_df: pd.DataFrame) -> pd.DataFrame:
    if len(meta_df) == 0:
        return pd.DataFrame()

    index = EvidenceIndex(evidence_df)

    rows = []
    for meta_dict in meta_df.to_dict("records"):
        subset = index.contributing(meta_dict)
        mapping_confidences = sorted(subset["mapping_confidence"].unique())
        quality_flags = sorted({f for flags in subset["quality_flags"] for f in _parse_flags(flags)})
        n_distinct_assays = index.n_distinct_assays(meta_dict)

        meta_dict["n_distinct_assays"] = n_distinct_assays
        components = compute_components(meta_dict, mapping_confidences, quality_flags)
        score = candidate_score(components)

        is_high_priority = (
            meta_dict["k_studies"] >= HIGH_PRIORITY_MIN_STUDIES
            and components["phenotype_relevance_score"] > 0.1
            and meta_dict["direction_consistency"] >= HIGH_PRIORITY_MIN_DIRECTION_CONSISTENCY
            and components["quality_score"] >= HIGH_PRIORITY_MIN_QUALITY_SCORE
            and significance_p_value(meta_dict) <= HIGH_PRIORITY_MAX_ADJUSTED_P
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
    return sort_candidates(out)


def sort_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    """Strongest tier first, then score descending within tier."""
    if len(candidates) == 0:
        return candidates
    tier_rank = candidates["tier"].map({t: i for i, t in enumerate(TIER_ORDER)}).fillna(len(TIER_ORDER))
    return (
        candidates.assign(_tier_rank=tier_rank)
        .sort_values(["_tier_rank", "score"], ascending=[True, False])
        .drop(columns="_tier_rank")
    )
