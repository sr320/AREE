"""Orchestrate meta-analysis over the harmonized evidence table.

Groups evidence records by (feature_id_standardized, phenotype, feature_type)
— never across feature types, since a log2FoldChange and a methylation
percent-difference are not on a comparable scale — and pools effect sizes
within each group. Records with mapping_confidence == 'unresolved' are
excluded from pooling (no stable identity to group on) but remain visible in
the raw evidence table.
"""
from __future__ import annotations

import pandas as pd

from common import EVIDENCE_TABLE_PATH, REPORTS_DIR

from .effect_sizes import effective_standard_error
from .pooling import dersimonian_laird

META_ANALYSIS_DIR = REPORTS_DIR / "meta_analysis"

GROUP_KEYS = ["feature_id_standardized", "phenotype", "feature_type"]


def _direction_consistency(effect_sizes: pd.Series) -> float:
    signs = effect_sizes.apply(lambda v: 1 if v > 0 else (-1 if v < 0 else 0))
    nonzero = signs[signs != 0]
    if len(nonzero) == 0:
        return 0.0
    majority_count = max((nonzero == 1).sum(), (nonzero == -1).sum())
    return float(majority_count / len(nonzero))


def run_meta_analysis(phenotype: str | None = None, feature_type: str | None = None) -> pd.DataFrame:
    if not EVIDENCE_TABLE_PATH.exists():
        raise FileNotFoundError(
            f"No evidence table at {EVIDENCE_TABLE_PATH}. Run `aree harmonize` for each study first."
        )
    df = pd.read_csv(EVIDENCE_TABLE_PATH, sep="\t")
    df = df[df["mapping_confidence"] != "unresolved"]
    df = df[df["feature_id_standardized"].notna() & df["effect_size"].notna()]

    if phenotype:
        df = df[df["phenotype"] == phenotype]
    if feature_type:
        df = df[df["feature_type"] == feature_type]

    results = []
    for keys, group in df.groupby(GROUP_KEYS, dropna=False):
        feature_id, pheno, ftype = keys
        ses = [
            effective_standard_error(row.effect_size, row.standard_error, row.p_value)
            for row in group.itertuples()
        ]
        usable = [(e, se) for e, se in zip(group["effect_size"], ses) if se is not None and se > 0]
        if not usable:
            continue
        yi = [u[0] for u in usable]
        sei = [u[1] for u in usable]
        pooled = dersimonian_laird(yi, sei)

        results.append({
            "feature_id_standardized": feature_id,
            "phenotype": pheno,
            "feature_type": ftype,
            "k_studies": group["study_id"].nunique(),
            "studies": "|".join(sorted(group["study_id"].unique())),
            "n_evidence_records": len(group),
            "total_sample_size": int(group.drop_duplicates(["study_id", "comparison_id"])["sample_size"].sum()),
            "pooled_effect": pooled.pooled_effect,
            "pooled_se": pooled.pooled_se,
            "ci_lower": pooled.ci_lower,
            "ci_upper": pooled.ci_upper,
            "z": pooled.z,
            "p_value": pooled.p_value,
            "q_statistic": pooled.q_statistic,
            "i_squared": pooled.i_squared,
            "tau_squared": pooled.tau_squared,
            "direction_consistency": _direction_consistency(group["effect_size"]),
            "distinct_tissues": group["tissue"].nunique(),
            "distinct_life_stages": group["life_stage"].nunique(),
            "distinct_stressors": "|".join(sorted(group["stressor"].unique())),
            "mapping_confidences": "|".join(sorted(group["mapping_confidence"].unique())),
            "quality_flags_union": "|".join(sorted({f for flags in group["quality_flags"] for f in _parse_flags(flags)})),
        })

    result_df = pd.DataFrame(results)
    if len(result_df):
        result_df = result_df.sort_values(["phenotype", "feature_type", "p_value"])
    return result_df


def _parse_flags(value) -> list[str]:
    if isinstance(value, list):
        return value
    if isinstance(value, str) and value.startswith("["):
        import ast
        return ast.literal_eval(value)
    return []


def write_meta_analysis(phenotype: str | None, feature_type: str | None) -> "pd.DataFrame | None":
    result = run_meta_analysis(phenotype, feature_type)
    META_ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    label = f"{phenotype or 'all_phenotypes'}_{feature_type or 'all_features'}"
    out_path = META_ANALYSIS_DIR / f"{label}_meta_analysis.tsv"
    result.to_csv(out_path, sep="\t", index=False)
    return result, out_path
