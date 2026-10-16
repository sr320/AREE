"""Cross-species evidence cards: one card per orthogroup, phenotype and feature type. WP6.

`aree build-cross-species-cards` runs the orthogroup-level meta-analysis, tiers
and scores every row (prioritize/cross_species.py), and writes, under
reports/evidence_cards/cross_species/:

* ``candidates.tsv``: every row, every tier, with its component breakdown and
  the gates it failed, so the full ranking is auditable;
* ``index.json``: the rows that got a card;
* one markdown card per row with a significant signal: the pooled
  cross-species adjusted p, or any one species' gene-level adjusted p, at or
  below the threshold. The second clause keeps discordant rows on the page:
  two species each significant in opposite directions pool to nothing.

Each card shows both (or all) species side by side: each species' own
gene-level estimate, then the pooled cross-species estimate and how much the
species disagree.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from common import load_vocab, species_term_for_taxid
from harmonize.core import load_evidence_table
from prioritize.cross_species import (
    TIER_ORDER,
    build_cross_species_candidates,
    species_estimates_for,
)
from prioritize.rank import _partition_value

from . import evidence_cards as _cards
from .evidence_cards import (
    EVIDENCE_CLASS_STATEMENT,
    NEXT_STEP_BY_EVIDENCE_CLASS,
    _fmt,
    _join,
    _slug,
)

CROSS_SPECIES_SUBDIR = "cross_species"

TIER_STATEMENT = {
    "cross_species_convergence": (
        "Every convergence gate passes: at least two species each replicate the association in two "
        "or more studies of their own, every species' estimate has the same sign, and the pooled "
        "cross-species estimate survives false-discovery control."
    ),
    "cross_species_emerging": (
        "The species agree in sign and the pooled estimate is significant, but at least one "
        "convergence gate fails (listed below). Treat it as a lead that needs replication, not as a "
        "conserved association."
    ),
    "cross_species_discordant": (
        "The species' estimates **disagree in sign**. This is shown, not dropped: it may be a real "
        "difference between species, or a difference between the studies (stressor, tissue, life "
        "stage, phenotype definition). The per-species table is the place to tell which."
    ),
    "cross_species_unsupported": (
        "The species agree in sign, but pooled across species the estimate is not significant."
    ),
}

NEXT_STEP_BY_TIER = {
    "cross_species_convergence": (
        "Test the association in an independent cohort of each species, then take the species with "
        "the strongest within-species evidence into functional follow-up."
    ),
    "cross_species_emerging": (
        "Replicate within each species that rests on a single study before treating the association "
        "as conserved."
    ),
    "cross_species_discordant": (
        "Compare the species' study designs before reading the reversal as biology. A study of both "
        "species under one matched design is what would resolve it."
    ),
    "cross_species_unsupported": (
        "No validation step is warranted from this evidence."
    ),
}

GATE_EXPLANATION = {
    "replication": "fewer than two species have k >= 2 studies of their own",
    "direction": "the species' estimates do not all share one sign",
    "significance": "the pooled cross-species adjusted p is above 0.05",
    "phenotype_relevance": "the weakest contributing evidence class is exposure-only",
    "quality": "quality flags across contributing records exceed the limit",
}


@dataclass
class CrossSpeciesCardResult:
    index: list
    n_candidates: int
    candidates_path: Path | None
    summary: dict


def cards_dir() -> Path:
    """Resolved at call time so tests can redirect EVIDENCE_CARDS_DIR."""
    return _cards.EVIDENCE_CARDS_DIR / CROSS_SPECIES_SUBDIR


def _species_label(taxid) -> str:
    term = species_term_for_taxid(taxid)
    return f"*{term['scientific_name']}*" if term else f"taxid {taxid}"


def _evidence_context(evidence: pd.DataFrame, ids: list[str]) -> dict:
    subset = evidence[evidence["evidence_id"].astype(str).isin(set(ids))]
    return {col: sorted({str(v) for v in subset[col].dropna()}) if col in subset.columns else []
            for col in ("feature_id_original", "tissue", "life_stage", "stressor", "study_id")}


def species_rows_for_card(candidate: dict, species_df: pd.DataFrame, evidence: pd.DataFrame) -> list[dict]:
    """One dict per species, with the evidence context a reader needs to compare them."""
    rows = []
    for r in species_estimates_for(species_df, candidate).to_dict("records"):
        ids = [i for i in str(r.get("contributing_evidence_ids") or "").split("|") if i]
        ctx = _evidence_context(evidence, ids)
        effect, se = float(r["pooled_effect"]), float(r["pooled_se"])
        rows.append({
            "species_taxid": r["species_taxid"], "species": _species_label(r["species_taxid"]),
            "gene": r["feature_id_standardized"], "original_ids": ctx["feature_id_original"],
            "k_studies": int(r["k_studies"]), "studies": str(r["studies"]),
            "tissues": ctx["tissue"], "life_stages": ctx["life_stage"], "stressors": ctx["stressor"],
            "sample_size": r.get("total_sample_size"),
            "effect": effect, "se": se, "ci_lower": effect - 1.959964 * se, "ci_upper": effect + 1.959964 * se,
            "adjusted_p_value": r.get("adjusted_p_value"),
        })
    return rows


def has_significant_signal(candidate: dict, species_rows: list[dict], max_adjusted_p: float) -> bool:
    q = candidate.get("adjusted_p_value")
    if q is not None and q == q and q <= max_adjusted_p:
        return True
    return any(r["adjusted_p_value"] is not None and r["adjusted_p_value"] == r["adjusted_p_value"]
               and r["adjusted_p_value"] <= max_adjusted_p for r in species_rows)


def _orthology_line(candidate: dict, run_ref: dict | None) -> str:
    if candidate["copy_class"] == "demo_crosswalk":
        return ("Orthology source: the demo crosswalk's illustrative orthogroup, one gene per species by "
                "construction. It is not from an OrthoFinder run.")
    run = run_ref["run"] if run_ref else "unknown (no active run)"
    copy = {
        "single_copy_all": "single-copy in every species in the run",
        "single_copy": "single-copy in each species that has it, absent from at least one",
        "multi_copy": ("**multi-copy**: at least one species has paralogs. Pooled only because each "
                       "species has evidence for exactly one gene; that gene may not be the true ortholog"),
    }.get(candidate["copy_class"], candidate["copy_class"])
    return f"Orthology source: OrthoFinder run `{run}`. Copy class: {copy}."


def _count(value) -> str:
    number = pd.to_numeric(value, errors="coerce")
    return "n/a" if pd.isna(number) or number <= 0 else str(int(number))


def _list(values: list) -> str:
    return ", ".join(values) if values else "not recorded"


def card_markdown(candidate: dict, species_rows: list[dict], phenotype_label: str, run_ref: dict | None) -> str:
    tier = candidate["tier"]
    evidence_class = candidate.get("evidence_class", "exposure_only")
    simulated = _partition_value(candidate["simulated"]) == "true"
    failed = [g for g in str(candidate.get("gates_failed") or "").split("|") if g]
    L = []
    L.append(f"# Cross-species evidence card: {candidate['orthogroup_id']} — {phenotype_label}")
    L.append("")
    L.append(f"**Tier:** `{tier}`  **Score:** {candidate['score']} / 100")
    L.append(f"**Evidence class:** `{evidence_class}`  **Copy class:** `{candidate['copy_class']}`  "
             f"**Species:** {candidate['n_species']}")
    L.append("")
    L.append("> A cross-species estimate says an association recurs in orthologous genes of different "
             "species. It is not a validated biomarker, and it does not show the genes do the same thing "
             "in each species. See Limitations below.")
    L.append("")
    L.append("## What this evidence does and does not show")
    L.append(TIER_STATEMENT[tier])
    L.append("")
    L.append(EVIDENCE_CLASS_STATEMENT.get(evidence_class, EVIDENCE_CLASS_STATEMENT["exposure_only"]))
    L.append("")
    L.append("## Orthology")
    L.append(f"- Orthogroup: `{candidate['orthogroup_id']}`")
    L.append(f"- {_orthology_line(candidate, run_ref)}")
    for r in species_rows:
        L.append(f"- {r['species']} (taxid `{r['species_taxid']}`): gene `{r['gene']}`, "
                 f"original ID(s) {_list(r['original_ids'])}")
    L.append(f"- Mapping confidence across contributing records: {candidate.get('mapping_confidences') or 'n/a'}")
    L.append(f"- Data origin: `{'simulated' if simulated else 'real'}`")
    L.append("")
    L.append("## Per-species estimates")
    L.append("Each species' estimate is its own gene-level random-effects pool, exactly as on its "
             "single-species card.")
    L.append("")
    L.append("| species | gene | studies (k) | tissue | life stage | stressor | n | effect | 95% CI | gene-level adj. p |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for r in species_rows:
        L.append(f"| {r['species']} | {r['gene']} | {r['studies'].replace('|', ', ')} ({r['k_studies']}) | "
                 f"{_list(r['tissues'])} | {_list(r['life_stages'])} | {_list(r['stressors'])} | "
                 f"{_count(r['sample_size'])} | "
                 f"{r['effect']:.3f} | {r['ci_lower']:.3f} to {r['ci_upper']:.3f} | "
                 f"{_fmt(r['adjusted_p_value'])} |")
    L.append("")
    L.append("## Cross-species pooled estimate")
    L.append(f"- Pooled effect across species: {candidate['pooled_effect']:.3f} "
             f"(95% CI {candidate['ci_lower']:.3f} to {candidate['ci_upper']:.3f})")
    L.append(f"- Pooled p-value: {_fmt(candidate['p_value'])}; BH-adjusted across the "
             f"{_fmt(candidate.get('n_tests_in_family'), 'd')} orthogroups tested for "
             f"`{candidate['phenotype']}` / `{candidate['feature_type']}`: {_fmt(candidate.get('adjusted_p_value'))}")
    L.append(f"- Heterogeneity between species: Q = {_fmt(candidate['q_between_species'])}, "
             f"I² = {_fmt(candidate['i_squared_between_species'], '.1f')}%, "
             f"τ² = {_fmt(candidate['tau_squared_between_species'])}")
    L.append(f"- Direction consistency across species: {candidate['direction_consistency_species']:.2f}")
    L.append(f"- Studies per species: {candidate['k_studies_by_species']}. Replication gate "
             f"(two or more species each with k ≥ 2): "
             f"{'met' if bool(candidate['meets_replication_gate']) else '**not met**'}")
    if failed:
        L.append(f"- Convergence gates failed: {_join([f'{g} ({GATE_EXPLANATION[g]})' for g in failed])}")
    if tier == "cross_species_discordant":
        L.append("- ⚠️ **Species disagree in direction.** Compare each species' confidence interval: an "
                 "interval spanning zero is noise around no effect, not a reversal.")
    L.append("")
    L.append("## Limitations")
    L.append("- Association, not mechanism, in every species shown.")
    if candidate["copy_class"] == "demo_crosswalk":
        L.append("- The orthogroup is illustrative, assigned by the demo crosswalk rather than inferred.")
    else:
        L.append("- Orthology is inferred from sequence. Shared ancestry does not guarantee shared function "
                 "or regulation.")
    designs = {k: {v for r in species_rows for v in r[k]} for k in ("tissues", "life_stages", "stressors")}
    differ = [label for key, label in (("stressors", "stressor"), ("tissues", "tissue"),
                                       ("life_stages", "life stage")) if len(designs[key]) > 1]
    if differ:
        L.append(f"- The species' studies differ in {_join(differ)}, so species and design are confounded "
                 "in the between-species comparison.")
    L.append(f"- Quality flags across contributing records: {candidate.get('quality_flags_union') or 'none recorded'}")
    if evidence_class != "resilience_associated":
        L.append(f"- At least one contributing study did not measure a resilience outcome, so the evidence "
                 f"class is `{evidence_class}`. A strong tier here is not a resilience claim.")
    if simulated:
        L.append("- This card contains simulated demo evidence; do not cite it as real evidence.")
    L.append("")
    L.append("## Recommended next validation step")
    phenotype_gap = NEXT_STEP_BY_EVIDENCE_CLASS.get(evidence_class)
    if phenotype_gap and tier in ("cross_species_convergence", "cross_species_emerging"):
        L.append(phenotype_gap)
        L.append("")
        L.append(f"Once a resilience outcome is available, the tier-appropriate step is: {NEXT_STEP_BY_TIER[tier]}")
    else:
        L.append(NEXT_STEP_BY_TIER[tier])
    L.append("")
    return "\n".join(L)


def card_slug(candidate: dict) -> str:
    origin = "simulated" if _partition_value(candidate["simulated"]) == "true" else "real"
    return _slug(f"{candidate['orthogroup_id']}_{candidate['phenotype']}_{candidate['feature_type']}_{origin}")


def build_cross_species_cards_report(
    phenotype: str | None = None,
    feature_type: str | None = None,
    max_adjusted_p: float = _cards.DEFAULT_MAX_ADJUSTED_P,
    all_cards: bool = False,
    include_multi_copy: bool = False,
) -> CrossSpeciesCardResult:
    from mappings.orthogroups import active_run_ref
    from meta_analysis.run import write_orthogroup_meta_analysis

    og_results, species_df, summary, _ = write_orthogroup_meta_analysis(
        phenotype, feature_type, include_multi_copy=include_multi_copy)
    evidence = load_evidence_table()
    candidates = build_cross_species_candidates(og_results, species_df, evidence)
    if len(candidates) == 0:
        return CrossSpeciesCardResult(index=[], n_candidates=0, candidates_path=None, summary=summary)

    out_dir = cards_dir()
    out_dir.mkdir(parents=True, exist_ok=True)
    phenotype_vocab = load_vocab("phenotype_ontology")
    run_ref = active_run_ref()

    card_files, written = [], []
    for candidate in candidates.to_dict("records"):
        species_rows = species_rows_for_card(candidate, species_df, evidence)
        if not all_cards and not has_significant_signal(candidate, species_rows, max_adjusted_p):
            card_files.append("")
            continue
        label = phenotype_vocab.get(candidate["phenotype"], {}).get("label", candidate["phenotype"])
        path = out_dir / f"{card_slug(candidate)}.md"
        path.write_text(card_markdown(candidate, species_rows, label, run_ref))
        card_files.append(str(path))
        written.append({k: candidate.get(k) for k in (
            "orthogroup_id", "phenotype", "feature_type", "simulated", "copy_class", "species_taxids",
            "tier", "evidence_class", "score", "k_studies_by_species", "adjusted_p_value", "gates_failed")}
            | {"card_file": str(path)})

    candidates = candidates.assign(card_file=card_files)
    candidates_path = out_dir / "candidates.tsv"
    candidates.to_csv(candidates_path, sep="\t", index=False)
    (out_dir / "index.json").write_text(json.dumps(written, indent=2, default=str))
    return CrossSpeciesCardResult(index=written, n_candidates=len(candidates),
                                  candidates_path=candidates_path, summary=summary)


__all__ = ["TIER_ORDER", "build_cross_species_cards_report", "card_markdown", "cards_dir"]
