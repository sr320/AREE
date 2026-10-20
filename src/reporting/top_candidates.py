"""Top-N candidate summary per phenotype, read from an already-ranked candidates.tsv.

`candidates.tsv` (written by `build_evidence_cards_report`) is one row per
feature x phenotype x feature_type x origin x species, globally sorted by tier
then score. This module re-groups that same ranking by phenotype so a reader
can see, per phenotype, which candidates lead the pack — without re-running
meta-analysis.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from common import REPORTS_DIR, load_vocab, species_term_for_taxid
from prioritize.rank import DEFAULT_CONTRAST_TYPE, TIER_ORDER

DEFAULT_CANDIDATES_PATH = REPORTS_DIR / "evidence_cards" / "candidates.tsv"
DEFAULT_OUT_PATH = REPORTS_DIR / "top_candidates_summary.md"

SUMMARY_COLUMNS = [
    "rank", "feature_id_standardized", "species", "feature_type", "tier", "evidence_class",
    "contrast_type", "score", "pooling_status", "k_studies", "total_sample_size", "direction_consistency",
    "adjusted_p_value", "context_replication", "n_supporting_layers",
    "mapping_confidences", "simulated", "card_file",
]


def load_candidates(path: Path = DEFAULT_CANDIDATES_PATH) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run `aree build-evidence-cards` first to rank candidates."
        )
    return pd.read_csv(path, sep="\t")


def _with_contrast_columns(candidates: pd.DataFrame) -> pd.DataFrame:
    """Candidates ranked before contrast types existed rest on exposure responses."""
    out = candidates.copy()
    if "contrast_type" not in out.columns:
        out["contrast_type"] = DEFAULT_CONTRAST_TYPE
    out["contrast_type"] = out["contrast_type"].fillna(DEFAULT_CONTRAST_TYPE)
    if "is_predictive_baseline" not in out.columns:
        out["is_predictive_baseline"] = False
    out["is_predictive_baseline"] = out["is_predictive_baseline"].astype(str).str.lower().eq("true")
    return out


def predictive_baseline_candidates(candidates: pd.DataFrame, n: int = 10) -> pd.DataFrame:
    """Top-n candidates resting on a baseline difference between resilient and
    susceptible groups, ranked by tier then score."""
    candidates = _with_contrast_columns(candidates)
    baseline = candidates[candidates["is_predictive_baseline"]]
    return _ranked(baseline).head(n).reset_index(drop=True)


def _ranked(candidates: pd.DataFrame) -> pd.DataFrame:
    tier_rank = candidates["tier"].map({t: i for i, t in enumerate(TIER_ORDER)}).fillna(len(TIER_ORDER))
    return (
        candidates.assign(_tier_rank=tier_rank)
        .sort_values(["_tier_rank", "score"], ascending=[True, False], na_position="last")
        .drop(columns="_tier_rank")
    )


def top_candidates_by_phenotype(candidates: pd.DataFrame, n: int = 10) -> dict:
    """phenotype -> top-n rows, ranked by tier then score (same order as candidates.tsv)."""
    tier_rank = candidates["tier"].map({t: i for i, t in enumerate(TIER_ORDER)}).fillna(len(TIER_ORDER))
    ranked = (
        candidates.assign(_tier_rank=tier_rank)
        .sort_values(["_tier_rank", "score"], ascending=[True, False])
        .drop(columns="_tier_rank")
    )
    return {
        phenotype: group.head(n).reset_index(drop=True)
        for phenotype, group in ranked.groupby("phenotype", sort=False)
    }


def _species_name(taxid) -> str:
    term = species_term_for_taxid(taxid)
    return term["scientific_name"] if term else f"taxid {taxid}"


def _phenotype_label(phenotype: str, vocab: dict) -> str:
    return vocab.get(phenotype, {}).get("label", phenotype)


def build_top_candidates_summary(
    n: int = 10,
    candidates_path: Path = DEFAULT_CANDIDATES_PATH,
    out_path: Path = DEFAULT_OUT_PATH,
    include_simulated: bool = False,
) -> Path:
    """Write a Markdown top-N-per-phenotype summary and return its path.

    Candidates from simulated demo studies are left out unless include_simulated
    is set: they demonstrate the pipeline and are not findings.
    """
    candidates = _with_contrast_columns(load_candidates(candidates_path))
    is_simulated = candidates["simulated"].astype(str).str.lower().eq("true")
    n_simulated_excluded = 0 if include_simulated else int(is_simulated.sum())
    if not include_simulated:
        candidates = candidates[~is_simulated].reset_index(drop=True)
    by_phenotype = top_candidates_by_phenotype(candidates, n=n)
    phenotype_vocab = load_vocab("phenotype_ontology")

    try:
        candidates_display = str(candidates_path.relative_to(REPORTS_DIR.parent))
    except ValueError:
        candidates_display = str(candidates_path)

    total = len(candidates)
    n_simulated_phenotypes = sorted(
        p for p, g in by_phenotype.items() if g["simulated"].astype(str).str.lower().eq("true").any()
    )

    lines = [
        "# Top candidate biomarkers by phenotype",
        "",
        f"Top {n} candidates per phenotype from `{candidates_display}` "
        f"({total} total ranked candidates), sorted by tier "
        f"({', '.join(TIER_ORDER)}) then transparent candidate score.",
        "",
        "A high score or top rank is not validation. `tier` reflects hard gates on replication, "
        "direction consistency, and significance (see docs/interpreting_candidate_scores.md); "
        "`evidence_class` states what the evidence is evidence *of* — resilience-associated, "
        "stress-response, disease-associated, or exposure-only — and only "
        "`resilience_associated` candidates bear on resilience directly.",
        "",
        "Rows with `pooling_status = unpooled` have no score: their source reports neither a "
        "standard error nor an unadjusted p-value (typically a significance-filtered DEG list), so "
        "they are listed for visibility, sorted after every pooled candidate, and never pooled.",
        "",
    ]
    if n_simulated_excluded:
        lines += [
            f"{n_simulated_excluded} candidates from simulated demo studies are left out. "
            "Pass `--include-simulated` to show them.",
            "",
        ]
    if not len(candidates):
        lines += ["No candidates from real studies yet.", ""]
    if n_simulated_phenotypes:
        lines += [
            f"Phenotypes with at least one simulated/demo candidate in this table: "
            f"{', '.join(n_simulated_phenotypes)}. Simulated evidence is clearly flagged in the "
            "`simulated` column below; treat it as a pipeline demonstration, not a finding.",
            "",
        ]

    baseline = predictive_baseline_candidates(candidates, n=n)
    n_baseline = int(candidates["is_predictive_baseline"].sum())
    lines += [
        "## Predictive baseline candidates",
        "",
        "Features that differ between groups of measured, differing resilience (line, family, "
        "population, ploidy or age) sampled under **common unstressed conditions** "
        "(`contrast_type = baseline_phenotype_difference`). This is the design for a marker that "
        "could be read before any stress; a response to stress is not. Tier still says how well "
        "replicated each one is: groups that differ in origin differ in more than resilience, so "
        "a single group pair is not enough (see docs/contrast_types.md).",
        "",
    ]
    if len(baseline):
        lines += [f"{len(baseline)} of {n_baseline} shown, all phenotypes.", "",
                  _summary_table(baseline, with_phenotype=True), ""]
    else:
        lines += ["No candidate rests on a baseline contrast with a resilience outcome yet.", ""]

    for phenotype, group in by_phenotype.items():
        label = _phenotype_label(phenotype, phenotype_vocab)
        lines.append(f"## {label} (`{phenotype}`)")
        lines.append("")
        lines.append(f"{len(group)} of {n} requested shown; "
                      f"{(candidates['phenotype'] == phenotype).sum()} total candidates for this phenotype.")
        lines.append("")
        lines.append(_summary_table(group))
        lines.append("")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")
    return out_path


def _summary_table(group: pd.DataFrame, with_phenotype: bool = False) -> str:
    table = group.copy()
    table.insert(0, "rank", range(1, len(table) + 1))
    table["card_file"] = table["card_file"].fillna("").apply(
        lambda p: Path(p).name if p else "(no significant signal — not rendered)"
    )
    # One phenotype can hold candidates from several species, and the same
    # feature ID may appear once per species; name the species on every row.
    table["species"] = table["species_taxid"].map(_species_name)
    # Candidates ranked before pooling_status existed are all pooled.
    if "pooling_status" not in table.columns:
        table["pooling_status"] = "pooled"
    table["pooling_status"] = table["pooling_status"].fillna("pooled")
    columns = list(SUMMARY_COLUMNS)
    if with_phenotype:
        columns.insert(columns.index("species"), "phenotype")
    return table[columns].to_markdown(index=False, floatfmt=".3g")
