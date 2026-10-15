"""Orthogroup-level (cross-species) meta-analysis. WP6 of docs/mussel_clam_integration_plan.md.

Gene-level meta-analysis (run.py) never pools across species: its group key
includes species_taxid. This module combines species through orthogroups, in
two stages:

1. **Within species.** For each orthogroup, phenotype and feature type, each
   species' estimate is that species' gene-level pooled estimate, exactly as
   run.py computed it, so every rule there (prespecified primary comparisons,
   duplicate-mapping exclusion, real standard errors only) still applies.
2. **Across species.** The per-species estimates are pooled with the same
   DerSimonian-Laird random-effects model. Its Q, I-squared and tau-squared
   now measure heterogeneity *between species*, which is reported, not hidden.

What may pool, and why:

* **Single-copy orthogroups only, by default.** In a single-copy orthogroup a
  species has at most one gene, so its estimate is one gene's estimate. In a
  multi-copy orthogroup, a species' "estimate" would be one paralog chosen from
  many, compared with another species' different paralog. Those orthogroups
  are counted and set aside unless `include_multi_copy` is given, and even then
  only when each species has evidence for exactly one gene.
* **Evidence from the active run only.** An orthogroup ID that the active run
  does not contain (evidence harmonized under an older run) is set aside.
* **The demo** carries its synthetic crosswalk's illustrative orthogroups, one
  gene each, classed ``demo_crosswalk``; simulated and real evidence never pool
  together, as at gene level.

Each result row also says whether it clears the WP6 replication gate: at least
two species, each with k >= 2 studies of its own. A cross-species estimate
resting on one study per species is reported but flagged.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .pooling import benjamini_hochberg, dersimonian_laird
from .run import _direction_consistency

POOLABLE_COPY_CLASSES = {"single_copy_all", "single_copy", "demo_crosswalk"}
FAMILY_KEYS = ["phenotype", "feature_type", "simulated"]
OG_GROUP_KEYS = ["orthogroup_id", "phenotype", "feature_type", "simulated"]

RESULT_COLUMNS = [
    "orthogroup_id", "copy_class", "phenotype", "feature_type", "simulated",
    "n_species", "species_taxids", "genes", "k_studies_by_species", "k_studies_total",
    "species_with_k2", "meets_replication_gate",
    "pooled_effect", "pooled_se", "ci_lower", "ci_upper", "z", "p_value",
    "adjusted_p_value", "n_tests_in_family",
    "q_between_species", "i_squared_between_species", "tau_squared_between_species",
    "direction_consistency_species", "species_effects",
]
SPECIES_COLUMNS = [
    "orthogroup_id", "phenotype", "feature_type", "simulated", "species_taxid",
    "feature_id_standardized", "k_studies", "studies", "pooled_effect", "pooled_se",
    "p_value", "adjusted_p_value",
]


def gene_to_orthogroup(evidence: pd.DataFrame) -> dict[tuple, str]:
    """(feature_id_standardized, simulated, species_taxid) -> orthogroup_id, from evidence.

    Every record of a gene carries the same orthogroup when it was harmonized
    under one run; a gene seen with two is a mix of runs and is an error.
    """
    ev = evidence[evidence["orthogroup_id"].notna() & evidence["feature_id_standardized"].notna()]
    keyed = ev.groupby(["feature_id_standardized", "simulated", "species_taxid"])["orthogroup_id"].unique()
    clashes = keyed[keyed.map(len) > 1]
    if len(clashes):
        raise ValueError(
            f"{len(clashes)} gene(s) carry two orthogroup IDs in the evidence table, e.g. "
            f"{clashes.index[0]}: {list(clashes.iloc[0])}. Re-harmonize under one orthogroup run.")
    return {k: v[0] for k, v in keyed.items()}


def _key(value):
    """Normalize simulated/species_taxid so gene-level results and evidence join."""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    text = str(value).strip()
    if text.lower() in ("true", "false"):
        return text.lower() == "true"
    try:
        return int(float(text))
    except ValueError:
        return text


def run_orthogroup_meta_analysis(
    gene_results: pd.DataFrame,
    evidence: pd.DataFrame,
    copy_classes: dict[str, str],
    include_multi_copy: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Pool gene-level estimates across species within orthogroups.

    `copy_classes` maps orthogroup_id -> copy class from the active run.
    Returns (cross-species results, per-species estimates, summary of exclusions).
    """
    summary = {"orthogroup_groups": 0, "single_species": 0, "multi_copy_excluded": 0,
               "paralogs_in_species_excluded": 0, "not_in_active_run_excluded": 0, "pooled": 0}
    if gene_results is None or len(gene_results) == 0:
        return pd.DataFrame(columns=RESULT_COLUMNS), pd.DataFrame(columns=SPECIES_COLUMNS), summary

    og_of = {(f, _key(s), _key(t)): og for (f, s, t), og in gene_to_orthogroup(evidence).items()}
    genes = gene_results.copy()
    genes["_sim"] = genes["simulated"].map(_key)
    genes["_taxid"] = genes["species_taxid"].map(_key)
    genes["orthogroup_id"] = [og_of.get((f, s, t)) for f, s, t in
                              zip(genes["feature_id_standardized"], genes["_sim"], genes["_taxid"])]
    genes = genes[genes["orthogroup_id"].notna()]

    rows, species_rows = [], []
    for (og, pheno, ftype, sim), g in genes.groupby(["orthogroup_id", "phenotype", "feature_type", "_sim"],
                                                   sort=True):
        summary["orthogroup_groups"] += 1
        if g["_taxid"].nunique() < 2:
            summary["single_species"] += 1
            continue
        if sim:
            # The demo crosswalk's illustrative IDs (OG0000001...) collide with
            # real OrthoFinder IDs, so simulated evidence never consults the run.
            copy_class = "demo_crosswalk"
        else:
            copy_class = copy_classes.get(og)
            if copy_class is None:
                summary["not_in_active_run_excluded"] += 1
                continue
        if copy_class not in POOLABLE_COPY_CLASSES and not include_multi_copy:
            summary["multi_copy_excluded"] += 1
            continue
        if (g.groupby("_taxid").size() > 1).any():
            summary["paralogs_in_species_excluded"] += 1
            continue

        g = g.sort_values("_taxid")
        effects = g["pooled_effect"].to_numpy(dtype=float)
        ses = g["pooled_se"].to_numpy(dtype=float)
        pooled = dersimonian_laird(effects, ses)
        k = g["k_studies"].astype(int).to_numpy()
        summary["pooled"] += 1
        rows.append({
            "orthogroup_id": og, "copy_class": copy_class, "phenotype": pheno,
            "feature_type": ftype, "simulated": sim,
            "n_species": len(g), "species_taxids": "|".join(str(t) for t in g["_taxid"]),
            "genes": "|".join(f"{t}:{f}" for t, f in zip(g["_taxid"], g["feature_id_standardized"])),
            "k_studies_by_species": "|".join(f"{t}:{n}" for t, n in zip(g["_taxid"], k)),
            "k_studies_total": int(k.sum()),
            "species_with_k2": int((k >= 2).sum()),
            "meets_replication_gate": bool((k >= 2).sum() >= 2),
            "pooled_effect": pooled.pooled_effect, "pooled_se": pooled.pooled_se,
            "ci_lower": pooled.ci_lower, "ci_upper": pooled.ci_upper,
            "z": pooled.z, "p_value": pooled.p_value,
            "q_between_species": pooled.q_statistic,
            "i_squared_between_species": pooled.i_squared,
            "tau_squared_between_species": pooled.tau_squared,
            # Share of species agreeing with the majority sign, as at gene level.
            "direction_consistency_species": _direction_consistency(effects),
            "species_effects": "|".join(f"{t}:{e:.3g}±{s:.2g}" for t, e, s in zip(g["_taxid"], effects, ses)),
        })
        for _, r in g.iterrows():
            species_rows.append({
                "orthogroup_id": og, "phenotype": pheno, "feature_type": ftype, "simulated": sim,
                "species_taxid": r["_taxid"], "feature_id_standardized": r["feature_id_standardized"],
                "k_studies": int(r["k_studies"]), "studies": r["studies"],
                "pooled_effect": r["pooled_effect"], "pooled_se": r["pooled_se"],
                "p_value": r["p_value"], "adjusted_p_value": r.get("adjusted_p_value"),
            })

    result = pd.DataFrame(rows)
    if len(result):
        adjusted = np.full(len(result), np.nan)
        n_tests = np.zeros(len(result), dtype=int)
        for positions in result.groupby(FAMILY_KEYS, sort=False).indices.values():
            adjusted[positions] = benjamini_hochberg(result["p_value"].to_numpy()[positions])
            n_tests[positions] = len(positions)
        result["adjusted_p_value"] = adjusted
        result["n_tests_in_family"] = n_tests
        result = result[RESULT_COLUMNS].sort_values(["phenotype", "feature_type", "p_value"])
    else:
        result = pd.DataFrame(columns=RESULT_COLUMNS)
    species = pd.DataFrame(species_rows, columns=SPECIES_COLUMNS)
    return result.reset_index(drop=True), species, summary
