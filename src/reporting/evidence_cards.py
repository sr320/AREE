"""Generate one evidence card per candidate biomarker.

Every registered candidate gets a card, regardless of tier — "emerging"
candidates are labeled as such rather than omitted, per docs/design.md's
requirement that every candidate (not just the strong ones) is auditable.
"""
from __future__ import annotations

import json
import re

import pandas as pd

from common import REPORTS_DIR, load_vocab
from harmonize.core import load_evidence_table
from meta_analysis.run import run_meta_analysis
from prioritize.rank import build_candidates

EVIDENCE_CARDS_DIR = REPORTS_DIR / "evidence_cards"

NEXT_STEP_BY_TIER = {
    "high_priority_cross_study": (
        "Prioritize for targeted validation (e.g. independent-cohort qPCR/Western confirmation) "
        "and functional follow-up."
    ),
    "multi_omics_convergence": (
        "Investigate the regulatory link across molecular layers (e.g. paired methylation/expression "
        "measurement in the same individuals) before committing to validation."
    ),
    "emerging": (
        "Requires replication in an independent study before further validation investment."
    ),
}


def _slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))


def _forest_points(evidence_df: pd.DataFrame, feature_id: str, phenotype: str, feature_type: str) -> list:
    subset = evidence_df[
        (evidence_df["feature_id_standardized"] == feature_id)
        & (evidence_df["phenotype"] == phenotype)
        & (evidence_df["feature_type"] == feature_type)
    ]
    points = []
    for _, r in subset.iterrows():
        points.append({
            "study_id": r["study_id"],
            "comparison_id": r["comparison_id"],
            "tissue": r["tissue"],
            "life_stage": r["life_stage"],
            "stressor": r["stressor"],
            "effect_size": r["effect_size"],
            "effect_size_type": r["effect_size_type"],
            "standard_error": r["standard_error"],
            "p_value": r["p_value"],
            "adjusted_p_value": r["adjusted_p_value"],
            "molecular_direction": r["molecular_direction"],
            "mapping_confidence": r["mapping_confidence"],
        })
    return points


def _other_omics_context(evidence_df: pd.DataFrame, feature_id: str, this_feature_type: str) -> list:
    """Evidence for this same standardized feature from OTHER assay/feature types,
    possibly under a different phenotype — this is the explicit, inspectable basis
    for the multi_omics_convergence tier, so a reader can see exactly which studies
    and phenotypes are being linked and judge for themselves whether that's a fair link."""
    other = evidence_df[
        (evidence_df["feature_id_standardized"] == feature_id)
        & (evidence_df["feature_type"] != this_feature_type)
    ]
    return other[[
        "study_id", "feature_type", "phenotype", "tissue", "molecular_direction", "effect_size", "mapping_confidence"
    ]].drop_duplicates().to_dict("records")


def _card_markdown(candidate: dict, forest_points: list, other_omics: list, phenotype_label: str) -> str:
    lines = []
    lines.append(f"# Evidence card: {candidate['feature_id_standardized']} — {phenotype_label}")
    lines.append("")
    lines.append(f"**Tier:** `{candidate['tier']}`  **Score:** {candidate['score']} / 100")
    lines.append("")
    lines.append("> AREE candidate scores are associations across available evidence, not validated biomarkers. "
                  "See Limitations below before acting on this card.")
    lines.append("")
    lines.append("## Candidate identifiers")
    lines.append(f"- Standardized feature ID: `{candidate['feature_id_standardized']}`")
    lines.append(f"- Feature type: `{candidate['feature_type']}`")
    lines.append(f"- Mapping confidence across contributing studies: {candidate['mapping_confidences']}")
    lines.append("")
    lines.append("## Evidence summary")
    lines.append(f"- Phenotype context: `{candidate['phenotype']}`")
    lines.append(f"- Independent studies supporting: {candidate['k_studies']} ({candidate['studies']})")
    lines.append(f"- Total biological sample size pooled: {candidate['total_sample_size']}")
    lines.append(f"- Pooled effect size: {candidate['pooled_effect']:.3f} "
                 f"(95% CI {candidate['ci_lower']:.3f} to {candidate['ci_upper']:.3f})")
    lines.append(f"- Pooled p-value: {candidate['p_value']:.3g}")
    lines.append(f"- Direction consistency across studies: {candidate['direction_consistency']:.2f}")
    lines.append(f"- Heterogeneity (I²): {candidate['i_squared']:.1f}%")
    if candidate["direction_consistency"] < 0.7 and candidate["k_studies"] >= 2:
        lines.append("- ⚠️ **Conflicting direction across studies** — treat as context-dependent, not a "
                      "stable biomarker, until the source of heterogeneity is resolved.")
    lines.append("")
    lines.append("## Per-study evidence (forest-plot data)")
    lines.append("| study | tissue | life stage | stressor | effect | SE | p | direction | mapping confidence |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for p in forest_points:
        se = f"{p['standard_error']:.3f}" if p["standard_error"] == p["standard_error"] and p["standard_error"] is not None else "n/a"
        pv = f"{p['p_value']:.3g}" if p["p_value"] == p["p_value"] and p["p_value"] is not None else "n/a"
        lines.append(
            f"| {p['study_id']} | {p['tissue']} | {p['life_stage']} | {p['stressor']} | "
            f"{p['effect_size']:.3f} | {se} | {pv} | {p['molecular_direction']} | {p['mapping_confidence']} |"
        )
    lines.append("")
    lines.append("## Multi-omics context")
    if other_omics:
        lines.append("Evidence for this same standardized feature from other assay types "
                      "(shown for transparency — these may reflect a different phenotype or tissue context; "
                      "review before treating as confirmatory):")
        lines.append("| study | assay/feature type | phenotype | tissue | direction | effect |")
        lines.append("|---|---|---|---|---|---|")
        for o in other_omics:
            lines.append(
                f"| {o['study_id']} | {o['feature_type']} | {o['phenotype']} | {o['tissue']} | "
                f"{o['molecular_direction']} | {o['effect_size']:.3f} |"
            )
    else:
        lines.append("No evidence from other assay/feature types for this standardized feature yet.")
    lines.append("")
    lines.append("## Limitations")
    lines.append(f"- Quality flags present across contributing evidence: {candidate['quality_flags_union'] or 'none recorded'}")
    lines.append("- This candidate reflects association, not confirmed mechanistic causation.")
    lines.append("- All contributing demo studies in this build use simulated data; do not cite this card as real evidence.")
    lines.append("")
    lines.append("## Recommended next validation step")
    lines.append(NEXT_STEP_BY_TIER.get(candidate["tier"], NEXT_STEP_BY_TIER["emerging"]))
    lines.append("")
    return "\n".join(lines)


def build_evidence_cards(phenotype: str | None = None, feature_type: str | None = None) -> list:
    evidence_df = load_evidence_table()
    if phenotype:
        evidence_df = evidence_df[evidence_df["phenotype"] == phenotype]

    feature_types = [feature_type] if feature_type else sorted(evidence_df["feature_type"].unique())

    all_candidates = []
    for ftype in feature_types:
        meta_df = run_meta_analysis(phenotype=phenotype, feature_type=ftype)
        if len(meta_df) == 0:
            continue
        candidates = build_candidates(meta_df, load_evidence_table())
        all_candidates.append(candidates)

    if not all_candidates:
        return []

    candidates_df = pd.concat(all_candidates, ignore_index=True)
    full_evidence = load_evidence_table()
    phenotype_vocab = load_vocab("phenotype_ontology")

    EVIDENCE_CARDS_DIR.mkdir(parents=True, exist_ok=True)
    index = []
    for _, candidate in candidates_df.iterrows():
        candidate = candidate.to_dict()
        forest_points = _forest_points(full_evidence, candidate["feature_id_standardized"], candidate["phenotype"], candidate["feature_type"])
        other_omics = _other_omics_context(full_evidence, candidate["feature_id_standardized"], candidate["feature_type"])
        phenotype_label = phenotype_vocab.get(candidate["phenotype"], {}).get("label", candidate["phenotype"])

        card_md = _card_markdown(candidate, forest_points, other_omics, phenotype_label)
        slug = _slug(f"{candidate['feature_id_standardized']}_{candidate['phenotype']}_{candidate['feature_type']}")
        card_path = EVIDENCE_CARDS_DIR / f"{slug}.md"
        card_path.write_text(card_md)

        index.append({
            "feature_id_standardized": candidate["feature_id_standardized"],
            "feature_type": candidate["feature_type"],
            "phenotype": candidate["phenotype"],
            "tier": candidate["tier"],
            "score": candidate["score"],
            "k_studies": candidate["k_studies"],
            "card_file": str(card_path),
        })

    (EVIDENCE_CARDS_DIR / "index.json").write_text(json.dumps(index, indent=2))
    return index
