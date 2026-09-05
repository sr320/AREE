"""Orchestrate meta-analysis over the harmonized evidence table.

Groups evidence records by standardized feature, phenotype, feature type,
simulation status, and species. Only statistically poolable records contribute
to the estimate or replication metadata. Records with unresolved identifiers
remain visible in the evidence table, while repeated within-study contrasts
fail closed because AREE does not yet model their covariance.

Every pooled p-value is then adjusted for multiple testing (Benjamini-Hochberg)
within its test family — the set of features pooled for one phenotype, feature
type, origin, and species. A genome-wide reanalysis contributes tens of
thousands of tests to a family, and an unadjusted pooled p is not evidence at
that scale.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm

from common import EVIDENCE_TABLE_PATH, REPORTS_DIR

from .effect_sizes import MIN_P
from .pooling import benjamini_hochberg, dersimonian_laird

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

# The multiple-testing family is everything in GROUP_KEYS except the feature
# itself. Because the CLI can only filter on phenotype and feature type — both
# family keys — the adjusted p-value of a given feature is the same whether the
# run covered one phenotype or all of them.
FAMILY_KEYS = ["phenotype", "feature_type", "simulated", "species_taxid"]

MAPPING_CONFIDENCE_RANK = {
    "exact": 5,
    "one_to_one_ortholog": 4,
    "inferred": 3,
    "many_to_one_ortholog": 2,
    "one_to_many_ortholog": 1,
    "unresolved": 0,
}

# Identifier columns must be read as text. Real NCBI GeneIDs are all digits, and
# pandas' default inference turned them into integers here while the evidence
# loader used by ranking and reporting kept them as strings — so no candidate
# row from a real study matched its own evidence records.
_ID_COLUMNS = {
    "evidence_id": str, "study_id": str, "comparison_id": str,
    "feature_id_original": str, "feature_id_standardized": str, "orthogroup_id": str,
}

RESULT_COLUMNS = [
    "feature_id_standardized", "phenotype", "feature_type", "simulated", "species_taxid",
    "k_studies", "studies", "n_evidence_records", "n_available_records",
    "n_excluded_unpoolable", "n_excluded_duplicate_mappings", "excluded_studies",
    "contributing_evidence_ids", "total_sample_size",
    "pooled_effect", "pooled_se", "ci_lower", "ci_upper", "z",
    "p_value", "adjusted_p_value", "n_tests_in_family",
    "q_statistic", "i_squared", "tau_squared", "direction_consistency",
    "distinct_tissues", "distinct_life_stages", "distinct_stressors",
    "mapping_confidences", "quality_flags_union",
]


def effective_standard_errors(evidence: pd.DataFrame) -> np.ndarray:
    """Vectorized twin of `effect_sizes.effective_standard_error`, for the whole table.

    Use a reported standard error when it is present and positive; otherwise back
    one out of the effect size and unadjusted p-value under the normal
    approximation (se = |effect| / Phi^-1(1 - p/2)). NaN where neither is
    possible. Computing this once for 50k rows is what removed the largest
    per-group cost from the meta-analysis.
    """
    effect = pd.to_numeric(evidence["effect_size"], errors="coerce").to_numpy(dtype=float)
    reported = pd.to_numeric(evidence["standard_error"], errors="coerce").to_numpy(dtype=float)
    p = pd.to_numeric(evidence["p_value"], errors="coerce").to_numpy(dtype=float)

    out = np.full(len(evidence), np.nan)
    use_reported = np.isfinite(reported) & (reported > 0)
    out[use_reported] = reported[use_reported]

    need = ~use_reported & np.isfinite(effect) & np.isfinite(p) & (effect != 0)
    if need.any():
        p_clipped = np.clip(p[need], MIN_P, 1.0 - 1e-16)
        z = norm.ppf(1 - p_clipped / 2)
        approx = np.where(z > 0, np.abs(effect[need]) / z, np.nan)
        out[need] = approx
    return out


def _describe_group(keys: tuple) -> str:
    feature_id, phenotype, feature_type, simulated, species_taxid = keys
    return (
        f"Group feature={feature_id!r}, phenotype={phenotype!r}, "
        f"feature_type={feature_type!r}, simulated={simulated!r}, "
        f"species_taxid={species_taxid!r}"
    )


def _deduplicate_identifier_collisions(
    study: np.ndarray, comparison: np.ndarray, rank: np.ndarray, original: np.ndarray, keys: tuple
) -> tuple[np.ndarray, int]:
    """Keep a uniquely best mapping when source IDs collapse within a comparison.

    A lower-confidence alias must not become a second biological effect. Equal-
    confidence collisions remain ambiguous and fail closed. Returns the positions
    (into the group arrays) to keep and how many rows were dropped.
    """
    if len(study) == 1:
        return np.array([0]), 0
    keep = []
    excluded = 0
    seen: dict = {}
    for pos in range(len(study)):
        seen.setdefault((study[pos], comparison[pos]), []).append(pos)
    for (study_id, comparison_id), positions in seen.items():
        if len(positions) == 1:
            keep.append(positions[0])
            continue
        ranks = rank[positions]
        best = [positions[i] for i in range(len(positions)) if ranks[i] == ranks.max()]
        if len(best) != 1:
            originals = ", ".join(sorted(str(original[i]) for i in positions))
            raise ValueError(
                "Multiple source identifiers with equal mapping confidence resolve to the same "
                f"meta-analysis feature. {_describe_group(keys)}; study={study_id!r}, "
                f"comparison={comparison_id!r}; source identifiers: {originals}. Resolve the "
                "identifier collision explicitly before pooling."
            )
        keep.append(best[0])
        excluded += len(positions) - 1
    return np.array(sorted(keep)), excluded


def _reject_correlated_within_study_effects(study: np.ndarray, comparison: np.ndarray, keys: tuple) -> None:
    """Fail closed when one study contributes multiple effects to one pool.

    AREE does not yet carry the covariance needed to combine contrasts that
    reuse samples or controls. Treating them as independent would understate
    uncertainty, so require curation/modeling rather than choosing a contrast
    or correlation silently.
    """
    if len(set(study)) == len(study):
        return
    by_study: dict = {}
    for study_id, comparison_id in zip(study, comparison):
        by_study.setdefault(study_id, set()).add(str(comparison_id))
    details = "; ".join(
        f"{study_id}: {', '.join(sorted(comparisons))}"
        for study_id, comparisons in sorted(by_study.items()) if len(comparisons) > 1
    )
    raise ValueError(
        "Meta-analysis cannot treat multiple comparisons from one study as independent. "
        f"{_describe_group(keys)} contains {details}. Select one prespecified "
        "comparison per study or implement a covariance-aware within-study model."
    )


def _direction_consistency(effects: np.ndarray) -> float:
    signs = np.sign(effects)
    nonzero = signs[signs != 0]
    if len(nonzero) == 0:
        return 0.0
    return float(max((nonzero > 0).sum(), (nonzero < 0).sum()) / len(nonzero))


def _load_evidence_for_pooling() -> pd.DataFrame:
    return pd.read_csv(EVIDENCE_TABLE_PATH, sep="\t", low_memory=False, dtype=_ID_COLUMNS)


def _joined(values) -> str:
    return "|".join(sorted({str(v) for v in values}))


def adjust_within_families(result_df: pd.DataFrame) -> pd.DataFrame:
    """Add BH-adjusted p-values and the family size to a meta-analysis table.

    Families are defined by FAMILY_KEYS. The adjustment is recomputed from the
    rows present, so it must be applied to a complete family — which
    `run_meta_analysis` guarantees, because its only filters are family keys.
    """
    result_df = result_df.copy()
    if len(result_df) == 0:
        result_df["adjusted_p_value"] = pd.Series(dtype=float)
        result_df["n_tests_in_family"] = pd.Series(dtype=int)
        return result_df
    adjusted = np.full(len(result_df), np.nan)
    n_tests = np.zeros(len(result_df), dtype=int)
    for positions in result_df.groupby(FAMILY_KEYS, dropna=False, sort=False).indices.values():
        adjusted[positions] = benjamini_hochberg(result_df["p_value"].to_numpy()[positions])
        n_tests[positions] = len(positions)
    result_df["adjusted_p_value"] = adjusted
    result_df["n_tests_in_family"] = n_tests
    return result_df


def run_meta_analysis(phenotype: str | None = None, feature_type: str | None = None) -> pd.DataFrame:
    if not EVIDENCE_TABLE_PATH.exists():
        raise FileNotFoundError(
            f"No evidence table at {EVIDENCE_TABLE_PATH}. Run `aree harmonize` for each study first."
        )
    df = _load_evidence_for_pooling()
    df = df[df["mapping_confidence"] != "unresolved"]
    df = df[df["feature_id_standardized"].notna() & df["effect_size"].notna()]

    if phenotype:
        df = df[df["phenotype"] == phenotype]
    if feature_type:
        df = df[df["feature_type"] == feature_type]
    df = df.reset_index(drop=True)
    if len(df) == 0:
        return pd.DataFrame()

    # Everything the per-group step needs, as flat arrays. Slicing a numpy array
    # by a handful of positions costs microseconds; slicing a DataFrame costs
    # milliseconds, and there are tens of thousands of groups.
    effect = pd.to_numeric(df["effect_size"], errors="coerce").to_numpy(dtype=float)
    se = effective_standard_errors(df)
    poolable_mask = np.isfinite(se) & (se > 0)
    sample_size = pd.to_numeric(df["sample_size"], errors="coerce").fillna(0).to_numpy()
    rank = df["mapping_confidence"].map(MAPPING_CONFIDENCE_RANK).fillna(-1).to_numpy()
    cols = {
        name: df[name].to_numpy(dtype=object)
        for name in ("study_id", "comparison_id", "evidence_id", "feature_id_original",
                     "tissue", "life_stage", "stressor", "mapping_confidence")
    }
    flag_cache: dict = {}
    flags = np.array(
        [tuple(flag_cache.setdefault(v, _parse_flags(v)) if isinstance(v, str) else _parse_flags(v))
         for v in df["quality_flags"].tolist()],
        dtype=object,
    )

    results = []
    for keys, positions in df.groupby(GROUP_KEYS, dropna=False, sort=True).indices.items():
        positions = np.asarray(positions)
        pool_pos = positions[poolable_mask[positions]]
        if len(pool_pos) == 0:
            continue
        unpool_pos = positions[~poolable_mask[positions]]

        keep, n_duplicate_mappings = _deduplicate_identifier_collisions(
            cols["study_id"][pool_pos], cols["comparison_id"][pool_pos],
            rank[pool_pos], cols["feature_id_original"][pool_pos], keys,
        )
        pool_pos = pool_pos[keep]
        _reject_correlated_within_study_effects(cols["study_id"][pool_pos], cols["comparison_id"][pool_pos], keys)

        pooled = dersimonian_laird(effect[pool_pos], se[pool_pos])
        feature_id, pheno, ftype, simulated, species_taxid = keys

        results.append({
            "feature_id_standardized": feature_id,
            "phenotype": pheno,
            "feature_type": ftype,
            "simulated": simulated,
            "species_taxid": species_taxid,
            "k_studies": len(set(cols["study_id"][pool_pos])),
            "studies": _joined(cols["study_id"][pool_pos]),
            "n_evidence_records": len(pool_pos),
            "n_available_records": len(positions),
            "n_excluded_unpoolable": len(unpool_pos),
            "n_excluded_duplicate_mappings": n_duplicate_mappings,
            "excluded_studies": _joined(cols["study_id"][unpool_pos]),
            "contributing_evidence_ids": _joined(cols["evidence_id"][pool_pos]),
            "total_sample_size": int(sample_size[pool_pos].sum()),
            "pooled_effect": pooled.pooled_effect,
            "pooled_se": pooled.pooled_se,
            "ci_lower": pooled.ci_lower,
            "ci_upper": pooled.ci_upper,
            "z": pooled.z,
            "p_value": pooled.p_value,
            "q_statistic": pooled.q_statistic,
            "i_squared": pooled.i_squared,
            "tau_squared": pooled.tau_squared,
            "direction_consistency": _direction_consistency(effect[pool_pos]),
            "distinct_tissues": len(set(cols["tissue"][pool_pos])),
            "distinct_life_stages": len(set(cols["life_stage"][pool_pos])),
            "distinct_stressors": _joined(cols["stressor"][pool_pos]),
            "mapping_confidences": _joined(cols["mapping_confidence"][pool_pos]),
            "quality_flags_union": _joined(f for fl in flags[pool_pos] for f in fl),
        })

    result_df = pd.DataFrame(results)
    if len(result_df):
        result_df = adjust_within_families(result_df)
        result_df = result_df[RESULT_COLUMNS]
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
