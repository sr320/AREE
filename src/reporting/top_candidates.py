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

from common import REPORTS_DIR, load_vocab
from prioritize.rank import TIER_ORDER

DEFAULT_CANDIDATES_PATH = REPORTS_DIR / "evidence_cards" / "candidates.tsv"
DEFAULT_OUT_PATH = REPORTS_DIR / "top_candidates_summary.md"

SUMMARY_COLUMNS = [
    "rank", "feature_id_standardized", "feature_type", "tier", "evidence_class",
    "score", "k_studies", "total_sample_size", "direction_consistency",
    "adjusted_p_value", "context_replication", "n_supporting_layers",
    "mapping_confidences", "simulated", "card_file",
]


def load_candidates(path: Path = DEFAULT_CANDIDATES_PATH) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run `aree build-evidence-cards` first to rank candidates."
        )
    return pd.read_csv(path, sep="\t")


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


def _phenotype_label(phenotype: str, vocab: dict) -> str:
    return vocab.get(phenotype, {}).get("label", phenotype)


def build_top_candidates_summary(
    n: int = 10,
    candidates_path: Path = DEFAULT_CANDIDATES_PATH,
    out_path: Path = DEFAULT_OUT_PATH,
) -> Path:
    """Write a Markdown top-N-per-phenotype summary and return its path."""
    candidates = load_candidates(candidates_path)
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
    ]
    if n_simulated_phenotypes:
        lines += [
            f"Phenotypes with at least one simulated/demo candidate in this table: "
            f"{', '.join(n_simulated_phenotypes)}. Simulated evidence is clearly flagged in the "
            "`simulated` column below; treat it as a pipeline demonstration, not a finding.",
            "",
        ]

    for phenotype, group in by_phenotype.items():
        label = _phenotype_label(phenotype, phenotype_vocab)
        lines.append(f"## {label} (`{phenotype}`)")
        lines.append("")
        lines.append(f"{len(group)} of {n} requested shown; "
                      f"{(candidates['phenotype'] == phenotype).sum()} total candidates for this phenotype.")
        lines.append("")
        table = group.copy()
        table.insert(0, "rank", range(1, len(table) + 1))
        table["card_file"] = table["card_file"].fillna("").apply(
            lambda p: Path(p).name if p else "(no significant signal — not rendered)"
        )
        table = table[SUMMARY_COLUMNS]
        lines.append(table.to_markdown(index=False, floatfmt=".3g"))
        lines.append("")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")
    return out_path
