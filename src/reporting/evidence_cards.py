"""Generate evidence cards for candidate biomarkers.

Every ranked candidate — including every "emerging" one — is written to
`reports/evidence_cards/candidates.tsv` with its full component breakdown, so
the complete ranking is auditable. A markdown card is rendered for the
candidates that carry at least one significant signal: a BH-adjusted pooled
p-value at or below the threshold, or a contributing study that itself
reported an adjusted p at or below it. The second clause is what keeps
conflicting-direction candidates on the page: two studies that each report a
significant effect in opposite directions pool to nothing, and that conflict
must be surfaced, not hidden.

A genome-wide pool produces a candidate per gene. Rendering ~30,000 cards for
genes that show no signal anywhere takes tens of minutes and buries the
reader; the threshold is what makes the card set readable. Pass
`all_cards=True` (CLI: `--all-cards`) to render every candidate regardless.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from common import REPORTS_DIR, load_vocab
from harmonize.core import load_evidence_table
from meta_analysis.run import run_meta_analysis
from prioritize.rank import (
    EvidenceIndex,
    _partition_value,
    add_partition_columns,
    build_candidates,
    sort_candidates,
)

EVIDENCE_CARDS_DIR = REPORTS_DIR / "evidence_cards"
DEFAULT_MAX_ADJUSTED_P = 0.05

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

FOREST_COLUMNS = [
    "study_id", "comparison_id", "tissue", "life_stage", "stressor", "effect_size",
    "effect_size_type", "standard_error", "p_value", "adjusted_p_value",
    "molecular_direction", "mapping_confidence",
]
OTHER_OMICS_COLUMNS = [
    "study_id", "feature_type", "_layer", "phenotype", "tissue", "molecular_direction",
    "effect_size", "adjusted_p_value", "mapping_confidence", "_significant",
]


def _slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))


def _fmt(value, spec: str = ".3g", missing: str = "n/a") -> str:
    if value is None or value != value:
        return missing
    return format(value, spec)


def _forest_points(index: EvidenceIndex, candidate: dict) -> list:
    subset = index.for_candidate(candidate)
    contributing_ids = set(str(candidate.get("contributing_evidence_ids") or "").split("|"))
    points = subset[FOREST_COLUMNS + ["evidence_id"]].to_dict("records")
    for point in points:
        point["included_in_pool"] = str(point.pop("evidence_id")) in contributing_ids
    return points


def _other_omics_context(index: EvidenceIndex, candidate: dict) -> list:
    """Evidence for this same standardized feature from OTHER molecular layers.

    Every row is shown, under any phenotype, so the reader can see all adjacent
    evidence. Each row is marked `counts_as_convergent_support` when it is for
    this candidate's own phenotype AND its study reported adjusted p <= 0.05 —
    those are the only rows the multi_omics_convergence gate credits."""
    other = index.other_layers(candidate)
    rows = other[OTHER_OMICS_COLUMNS].drop_duplicates().to_dict("records")
    for row in rows:
        row["layer"] = row.pop("_layer")
        row["significant"] = bool(row.pop("_significant"))
        row["same_phenotype"] = row["phenotype"] == candidate["phenotype"]
        row["counts_as_convergent_support"] = row["same_phenotype"] and row["significant"]
    return rows


def has_significant_signal(candidate: dict, forest_points: list, max_adjusted_p: float) -> bool:
    """Does anything about this candidate clear the significance threshold?

    Either the pooled estimate does (BH-adjusted within its family), or at
    least one contributing study reported an adjusted p that does.
    """
    pooled_q = candidate.get("adjusted_p_value")
    if pooled_q is not None and pooled_q == pooled_q and pooled_q <= max_adjusted_p:
        return True
    return any(
        p["adjusted_p_value"] is not None
        and p["adjusted_p_value"] == p["adjusted_p_value"]
        and p["adjusted_p_value"] <= max_adjusted_p
        for p in forest_points
    )


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
    lines.append(f"- Species taxid: `{candidate['species_taxid']}`")
    simulated = _partition_value(candidate["simulated"]) == "true"
    lines.append(f"- Data origin: `{'simulated' if simulated else 'real'}`")
    lines.append(f"- Mapping confidence across contributing studies: {candidate['mapping_confidences']}")
    lines.append("")
    lines.append("## Evidence summary")
    lines.append(f"- Phenotype context: `{candidate['phenotype']}`")
    lines.append(f"- Independent studies supporting: {candidate['k_studies']} ({candidate['studies']})")
    lines.append(f"- Total biological sample size pooled: {candidate['total_sample_size']}")
    lines.append(f"- Pooled effect size: {candidate['pooled_effect']:.3f} "
                 f"(95% CI {candidate['ci_lower']:.3f} to {candidate['ci_upper']:.3f})")
    lines.append(f"- Pooled p-value: {_fmt(candidate['p_value'])}")
    n_tests = candidate.get("n_tests_in_family")
    lines.append(
        f"- Pooled p-value, BH-adjusted across the {_fmt(n_tests, 'd')} features tested for "
        f"`{candidate['phenotype']}` / `{candidate['feature_type']}` in this partition: "
        f"{_fmt(candidate.get('adjusted_p_value'))}"
    )
    lines.append(f"- Direction consistency across studies: {candidate['direction_consistency']:.2f}")
    lines.append(f"- Heterogeneity (I²): {candidate['i_squared']:.1f}%")
    if candidate["direction_consistency"] < 0.7 and candidate["k_studies"] >= 2:
        lines.append("- ⚠️ **Conflicting direction across studies** — treat as context-dependent, not a "
                      "stable biomarker, until the source of heterogeneity is resolved.")
    lines.append("")
    lines.append("## Per-study evidence (forest-plot data)")
    lines.append("| study | tissue | life stage | stressor | effect | SE | p | adj. p | pooled | direction | mapping confidence |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for p in forest_points:
        lines.append(
            f"| {p['study_id']} | {p['tissue']} | {p['life_stage']} | {p['stressor']} | "
            f"{p['effect_size']:.3f} | {_fmt(p['standard_error'], '.3f')} | {_fmt(p['p_value'])} | "
            f"{_fmt(p['adjusted_p_value'])} | "
            f"{'yes' if p['included_in_pool'] else 'no'} | {p['molecular_direction']} | "
            f"{p['mapping_confidence']} |"
        )
    lines.append("")
    lines.append("## Multi-omics context")
    own_layer = candidate.get("molecular_layer", candidate["feature_type"])
    supporting = [x for x in str(candidate.get("supporting_layers") or "").split("|") if x]
    lines.append(
        f"This candidate's own layer is `{own_layer}`. Layers with a significant record "
        f"(adjusted p ≤ 0.05) for this feature under `{candidate['phenotype']}`: "
        f"{', '.join(f'`{x}`' for x in supporting) if supporting else 'none'}."
    )
    if candidate.get("is_multi_omics_convergence"):
        lines.append(
            "**How the layers were linked:** each layer's records were mapped to the same "
            f"standardized identifier `{candidate['feature_id_standardized']}` "
            f"(mapping confidence: {candidate['mapping_confidences']}); the link is shared gene "
            "identity, not a measured regulatory relationship. Direction across layers is shown "
            "below and is not required to agree."
        )
    if other_omics:
        lines.append("")
        lines.append("Evidence for this same standardized feature from other molecular layers, under any "
                      "phenotype. Only rows marked **yes** in the last column count toward convergence; the "
                      "rest are shown so adjacent evidence is not hidden:")
        lines.append("| study | layer | feature type | phenotype | tissue | direction | effect | adj. p | same phenotype | significant | counts |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for o in other_omics:
            lines.append(
                f"| {o['study_id']} | {o['layer']} | {o['feature_type']} | {o['phenotype']} | {o['tissue']} | "
                f"{o['molecular_direction']} | {o['effect_size']:.3f} | {_fmt(o['adjusted_p_value'])} | "
                f"{'yes' if o['same_phenotype'] else 'no'} | {'yes' if o['significant'] else 'no'} | "
                f"{'**yes**' if o['counts_as_convergent_support'] else 'no'} |"
            )
    else:
        lines.append("")
        lines.append("No evidence from other molecular layers for this standardized feature yet.")
    lines.append("")
    lines.append("## Limitations")
    lines.append(f"- Quality flags present across contributing evidence: {candidate['quality_flags_union'] or 'none recorded'}")
    lines.append("- This candidate reflects association, not confirmed mechanistic causation.")
    if simulated:
        lines.append("- This card contains simulated demo evidence; do not cite it as real evidence.")
    lines.append("")
    lines.append("## Recommended next validation step")
    lines.append(NEXT_STEP_BY_TIER.get(candidate["tier"], NEXT_STEP_BY_TIER["emerging"]))
    lines.append("")
    return "\n".join(lines)


def _card_slug(candidate: dict) -> str:
    origin = "simulated" if _partition_value(candidate["simulated"]) == "true" else "real"
    return _slug(
        f"{candidate['feature_id_standardized']}_{candidate['phenotype']}_"
        f"{candidate['feature_type']}_{origin}_taxid{candidate['species_taxid']}"
    )


def rank_all_candidates(evidence_df: pd.DataFrame, phenotype: str | None, feature_type: str | None) -> pd.DataFrame:
    """Meta-analyze and rank every candidate matching the filters, one evidence pass."""
    scope = evidence_df
    if phenotype:
        scope = scope[scope["phenotype"] == phenotype]
    feature_types = [feature_type] if feature_type else sorted(scope["feature_type"].dropna().unique())

    frames = []
    for ftype in feature_types:
        meta_df = run_meta_analysis(phenotype=phenotype, feature_type=ftype)
        if len(meta_df):
            frames.append(build_candidates(meta_df, evidence_df))
    if not frames:
        return pd.DataFrame()
    return sort_candidates(pd.concat(frames, ignore_index=True))


@dataclass
class CardBuildResult:
    index: list                  # one entry per card written (also saved as index.json)
    n_candidates: int            # every ranked candidate, written or not
    candidates_path: Path | None # reports/evidence_cards/candidates.tsv, or None if no candidates


def build_evidence_cards(
    phenotype: str | None = None,
    feature_type: str | None = None,
    max_adjusted_p: float = DEFAULT_MAX_ADJUSTED_P,
    all_cards: bool = False,
) -> list:
    """Rank candidates, write `candidates.tsv` for all of them, and render a
    markdown card for each one that carries a significant signal.

    Returns the index of cards written (also saved as `index.json`).
    """
    return build_evidence_cards_report(phenotype, feature_type, max_adjusted_p, all_cards).index


def build_evidence_cards_report(
    phenotype: str | None = None,
    feature_type: str | None = None,
    max_adjusted_p: float = DEFAULT_MAX_ADJUSTED_P,
    all_cards: bool = False,
) -> CardBuildResult:
    evidence = add_partition_columns(load_evidence_table())
    candidates_df = rank_all_candidates(evidence, phenotype, feature_type)
    if len(candidates_df) == 0:
        return CardBuildResult(index=[], n_candidates=0, candidates_path=None)

    EVIDENCE_CARDS_DIR.mkdir(parents=True, exist_ok=True)
    index = EvidenceIndex(evidence)
    phenotype_vocab = load_vocab("phenotype_ontology")

    card_files = []
    written = []
    for candidate in candidates_df.to_dict("records"):
        forest_points = _forest_points(index, candidate)
        if not all_cards and not has_significant_signal(candidate, forest_points, max_adjusted_p):
            card_files.append("")
            continue

        other_omics = _other_omics_context(index, candidate)
        phenotype_label = phenotype_vocab.get(candidate["phenotype"], {}).get("label", candidate["phenotype"])
        card_path = EVIDENCE_CARDS_DIR / f"{_card_slug(candidate)}.md"
        card_path.write_text(_card_markdown(candidate, forest_points, other_omics, phenotype_label))
        card_files.append(str(card_path))
        written.append({
            "feature_id_standardized": candidate["feature_id_standardized"],
            "feature_type": candidate["feature_type"],
            "phenotype": candidate["phenotype"],
            "simulated": candidate["simulated"],
            "species_taxid": candidate["species_taxid"],
            "tier": candidate["tier"],
            "score": candidate["score"],
            "k_studies": candidate["k_studies"],
            "adjusted_p_value": candidate.get("adjusted_p_value"),
            "card_file": str(card_path),
        })

    candidates_df = candidates_df.assign(card_file=card_files)
    candidates_path = EVIDENCE_CARDS_DIR / "candidates.tsv"
    candidates_df.to_csv(candidates_path, sep="\t", index=False)
    (EVIDENCE_CARDS_DIR / "index.json").write_text(json.dumps(written, indent=2, default=str))
    return CardBuildResult(index=written, n_candidates=len(candidates_df), candidates_path=candidates_path)
