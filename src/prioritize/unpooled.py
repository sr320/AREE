"""Candidates that rest only on evidence meta-analysis cannot pool.

Inverse-variance pooling needs a standard error, reported or backed out of an
unadjusted p-value. Some public sources publish neither — typically an
author-curated, significance-filtered DEG list with log2FC and padj only
(HESSER2024_VCOR). Backing a standard error out of padj, or pooling a
significance-selected subset at all, would bias pooled effects upward and
inflate apparent replication, so those records stay out of every pooled
estimate (see run_meta_analysis).

Leaving them out of the pool must not make them disappear. A feature whose only
evidence for a phenotype is unpoolable would otherwise get no row in
candidates.tsv — and for HESSER2024_VCOR that silently dropped the only real
`disease_resistance` comparison in the registry. This module gives each such
feature an `emerging` row with `pooling_status = unpooled` and no score: the
transparent score is defined on pooled quantities, and a number computed from a
significance-filtered record would rank it above honestly pooled candidates.
These rows sort after every scored candidate.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from meta_analysis.run import _non_primary_comparisons, effective_standard_errors

from .rank import (
    CANDIDATE_GROUP_COLUMNS,
    add_partition_columns,
    candidate_key,
    conservative_relevance,
    evidence_class,
    layer_of,
)

POOLING_STATUS_POOLED = "pooled"
POOLING_STATUS_UNPOOLED = "unpooled"

UNPOOLED_REASON = "no_standard_error_or_unadjusted_p"


def _parse_flags(value) -> list:
    if isinstance(value, list):
        return value
    if isinstance(value, str) and value.startswith("["):
        import ast
        return ast.literal_eval(value)
    return []


def unpoolable_records(evidence_df: pd.DataFrame) -> pd.DataFrame:
    """Records eligible for pooling in every respect except a usable standard error.

    Mirrors run_meta_analysis's eligibility: a resolved identifier, an effect
    size, and a comparison the study has not prespecified out of pooling.
    """
    evidence = add_partition_columns(evidence_df)
    evidence = evidence[
        (evidence["mapping_confidence"] != "unresolved")
        & evidence["feature_id_standardized"].notna()
        & pd.to_numeric(evidence["effect_size"], errors="coerce").notna()
    ]
    if len(evidence) == 0:
        return evidence
    se = effective_standard_errors(evidence)
    unpoolable = ~(np.isfinite(se) & (se > 0))
    non_primary = _non_primary_comparisons()
    primary = np.array(
        [(s, c) not in non_primary for s, c in zip(evidence["study_id"], evidence["comparison_id"])],
        dtype=bool,
    )
    out = evidence[unpoolable & primary].copy()
    out["feature_id_standardized"] = out["feature_id_standardized"].astype(str)
    return out


def build_unpooled_candidates(evidence_df: pd.DataFrame, pooled_candidates: pd.DataFrame,
                              phenotype: str | None = None, feature_type: str | None = None) -> pd.DataFrame:
    """One `emerging`, unscored row per feature/phenotype/feature type/partition
    whose evidence is entirely unpoolable."""
    records = unpoolable_records(evidence_df)
    if phenotype:
        records = records[records["phenotype"] == phenotype]
    if feature_type:
        records = records[records["feature_type"] == feature_type]
    if len(records) == 0:
        return pd.DataFrame()

    pooled_keys = (
        {candidate_key(row) for row in pooled_candidates.to_dict("records")}
        if len(pooled_candidates) else set()
    )

    rows = []
    for _, group in records.groupby(CANDIDATE_GROUP_COLUMNS, sort=True):
        first = group.iloc[0]
        meta_like = {
            "feature_id_standardized": first["feature_id_standardized"],
            "phenotype": first["phenotype"],
            "feature_type": first["feature_type"],
            "simulated": first["simulated"],
            "species_taxid": first["species_taxid"],
        }
        if candidate_key(meta_like) in pooled_keys:
            continue
        classifications = (group["resilience_classification"]
                           if "resilience_classification" in group.columns else [])
        relevance = conservative_relevance(first["phenotype"], classifications)
        effects = pd.to_numeric(group["effect_size"], errors="coerce")
        signs = np.sign(effects[effects != 0])
        flags = sorted({f for v in group["quality_flags"] for f in _parse_flags(v)})
        studies = "|".join(sorted(set(group["study_id"].astype(str))))
        rows.append({
            **meta_like,
            "k_studies": 0,
            "studies": "",
            "n_evidence_records": 0,
            "n_available_records": len(group),
            "n_excluded_unpoolable": len(group),
            "excluded_studies": studies,
            "unpooled_evidence_ids": "|".join(sorted(group["evidence_id"].astype(str))),
            "unpooled_reason": UNPOOLED_REASON,
            "total_sample_size": int(pd.to_numeric(group["sample_size"], errors="coerce").fillna(0).sum()),
            "direction_consistency": (float(max((signs > 0).sum(), (signs < 0).sum()) / len(signs))
                                      if len(signs) else np.nan),
            "distinct_tissues": group["tissue"].nunique(),
            "distinct_life_stages": group["life_stage"].nunique(),
            "distinct_stressors": "|".join(sorted(set(group["stressor"].astype(str)))),
            "mapping_confidences": "|".join(sorted(set(group["mapping_confidence"].astype(str)))),
            "quality_flags_union": "|".join(flags),
            "score": np.nan,
            "tier": "emerging",
            "evidence_class": evidence_class(meta_like, relevance),
            "context_replication": "single_study",
            "pooling_status": POOLING_STATUS_UNPOOLED,
            "orthogroup_id": ("|".join(sorted(group["orthogroup_id"].dropna().astype(str).unique()))
                              if "orthogroup_id" in group.columns else ""),
            "molecular_layer": layer_of(first["feature_type"]),
        })
    return pd.DataFrame(rows)
