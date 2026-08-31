"""Orchestrate meta-analysis over the harmonized evidence table.

Groups evidence records by standardized feature, phenotype, feature type,
simulation status, and species. Only statistically poolable records contribute
to the estimate or replication metadata. Records with unresolved identifiers
remain visible in the evidence table, while repeated within-study contrasts
fail closed because AREE does not yet model their covariance.
"""
from __future__ import annotations

import pandas as pd

from common import EVIDENCE_TABLE_PATH, REPORTS_DIR

from .effect_sizes import effective_standard_error
from .pooling import dersimonian_laird

META_ANALYSIS_DIR = REPORTS_DIR / "meta_analysis"

# `simulated` is part of the grouping key so that a pooled estimate can never mix
# fabricated demo evidence with evidence from a real study.
#
# `species_taxid` is the canonical species identity rather than the reported name:
# grouping on the free-text name would both split one species across its accepted
# synonyms (Crassostrea gigas / Magallana gigas) and, once a second species is
# registered, silently pool across species. Cross-species evidence is only ever
# meant to combine through orthogroups, which are not yet populated.
GROUP_KEYS = ["feature_id_standardized", "phenotype", "feature_type", "simulated", "species_taxid"]
MAPPING_CONFIDENCE_RANK = {
    "exact": 5,
    "one_to_one_ortholog": 4,
    "inferred": 3,
    "many_to_one_ortholog": 2,
    "one_to_many_ortholog": 1,
    "unresolved": 0,
}


def _poolable_group(group: pd.DataFrame) -> pd.DataFrame:
    """Return only records that can actually enter inverse-variance pooling."""
    poolable = group.copy()
    poolable["_effective_se"] = [
        effective_standard_error(row.effect_size, row.standard_error, row.p_value)
        for row in group.itertuples()
    ]
    return poolable[
        poolable["_effective_se"].notna() & (poolable["_effective_se"] > 0)
    ].copy()


def _deduplicate_identifier_collisions(group: pd.DataFrame, keys: tuple) -> tuple[pd.DataFrame, int]:
    """Keep a uniquely best mapping when source IDs collapse within a comparison.

    A lower-confidence alias must not become a second biological effect. Equal-
    confidence collisions remain ambiguous and fail closed.
    """
    selected = []
    excluded_count = 0
    for (study_id, comparison_id), comparison_group in group.groupby(
        ["study_id", "comparison_id"], sort=False
    ):
        if len(comparison_group) == 1:
            selected.append(comparison_group.index[0])
            continue
        ranks = comparison_group["mapping_confidence"].map(MAPPING_CONFIDENCE_RANK).fillna(-1)
        best = comparison_group[ranks == ranks.max()]
        if len(best) != 1:
            feature_id, phenotype, feature_type, simulated, species_taxid = keys
            originals = ", ".join(sorted(comparison_group["feature_id_original"].astype(str)))
            raise ValueError(
                "Multiple source identifiers with equal mapping confidence resolve to the same "
                f"meta-analysis feature. Group feature={feature_id!r}, phenotype={phenotype!r}, "
                f"feature_type={feature_type!r}, simulated={simulated!r}, "
                f"species_taxid={species_taxid!r}; study={study_id!r}, "
                f"comparison={comparison_id!r}; source identifiers: {originals}. Resolve the "
                "identifier collision explicitly before pooling."
            )
        selected.append(best.index[0])
        excluded_count += len(comparison_group) - 1
    return group.loc[selected].copy(), excluded_count


def _reject_correlated_within_study_effects(group: pd.DataFrame, keys: tuple) -> None:
    """Fail closed when one study contributes multiple effects to one pool.

    AREE does not yet carry the covariance needed to combine contrasts that
    reuse samples or controls. Treating them as independent would understate
    uncertainty, so require curation/modeling rather than choosing a contrast
    or correlation silently.
    """
    duplicated = group[group.duplicated("study_id", keep=False)]
    if duplicated.empty:
        return
    details = "; ".join(
        f"{study_id}: {', '.join(sorted(study_group['comparison_id'].astype(str).unique()))}"
        for study_id, study_group in duplicated.groupby("study_id")
    )
    feature_id, phenotype, feature_type, simulated, species_taxid = keys
    raise ValueError(
        "Meta-analysis cannot treat multiple comparisons from one study as independent. "
        f"Group feature={feature_id!r}, phenotype={phenotype!r}, "
        f"feature_type={feature_type!r}, simulated={simulated!r}, "
        f"species_taxid={species_taxid!r} contains {details}. Select one prespecified "
        "comparison per study or implement a covariance-aware within-study model."
    )


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
        feature_id, pheno, ftype, simulated, species_taxid = keys
        poolable_candidates = _poolable_group(group)
        unpoolable = group.loc[~group.index.isin(poolable_candidates.index)]
        if poolable_candidates.empty:
            continue
        poolable, n_duplicate_mappings = _deduplicate_identifier_collisions(poolable_candidates, keys)
        _reject_correlated_within_study_effects(poolable, keys)

        yi = poolable["effect_size"].astype(float).tolist()
        sei = poolable["_effective_se"].astype(float).tolist()
        pooled = dersimonian_laird(yi, sei)

        results.append({
            "feature_id_standardized": feature_id,
            "phenotype": pheno,
            "feature_type": ftype,
            "simulated": simulated,
            "species_taxid": species_taxid,
            "k_studies": poolable["study_id"].nunique(),
            "studies": "|".join(sorted(poolable["study_id"].unique())),
            "n_evidence_records": len(poolable),
            "n_available_records": len(group),
            "n_excluded_unpoolable": len(unpoolable),
            "n_excluded_duplicate_mappings": n_duplicate_mappings,
            "excluded_studies": "|".join(sorted(unpoolable["study_id"].unique())),
            "contributing_evidence_ids": "|".join(sorted(poolable["evidence_id"].astype(str))),
            "total_sample_size": int(poolable["sample_size"].sum()),
            "pooled_effect": pooled.pooled_effect,
            "pooled_se": pooled.pooled_se,
            "ci_lower": pooled.ci_lower,
            "ci_upper": pooled.ci_upper,
            "z": pooled.z,
            "p_value": pooled.p_value,
            "q_statistic": pooled.q_statistic,
            "i_squared": pooled.i_squared,
            "tau_squared": pooled.tau_squared,
            "direction_consistency": _direction_consistency(poolable["effect_size"]),
            "distinct_tissues": poolable["tissue"].nunique(),
            "distinct_life_stages": poolable["life_stage"].nunique(),
            "distinct_stressors": "|".join(sorted(poolable["stressor"].unique())),
            "mapping_confidences": "|".join(sorted(poolable["mapping_confidence"].unique())),
            "quality_flags_union": "|".join(sorted({f for flags in poolable["quality_flags"] for f in _parse_flags(flags)})),
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


def write_meta_analysis(phenotype: str | None, feature_type: str | None) -> pd.DataFrame | None:
    result = run_meta_analysis(phenotype, feature_type)
    META_ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    label = f"{phenotype or 'all_phenotypes'}_{feature_type or 'all_features'}"
    out_path = META_ANALYSIS_DIR / f"{label}_meta_analysis.tsv"
    result.to_csv(out_path, sep="\t", index=False)
    return result, out_path
