"""Read-only data access for the Streamlit app.

The app never mutates the registry or evidence tables — it only reads the
artifacts produced by the `aree` CLI. If those artifacts are missing, the app
tells the user which command to run rather than failing silently.

Simulated demo studies and their evidence are left out unless a loader is
called with include_simulated=True (the sidebar toggle in main.py).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from common import (  # noqa: E402
    DEMO_REGISTRY_CSV,
    EVIDENCE_TABLE_PATH,
    STUDY_REGISTRY_CSV,
    load_vocab,
    load_yaml,
    resolve_species,
    species_term_for_taxid,
    study_files,
)
from harmonize.core import load_evidence_table  # noqa: E402
from meta_analysis.run import run_meta_analysis  # noqa: E402
from prioritize.rank import build_candidates  # noqa: E402
from reporting.species import species_summary  # noqa: E402


def registry_exists() -> bool:
    return STUDY_REGISTRY_CSV.exists() and STUDY_REGISTRY_CSV.stat().st_size > 0


def evidence_exists() -> bool:
    return EVIDENCE_TABLE_PATH.exists()


def _is_simulated(df: pd.DataFrame) -> pd.Series:
    return df["simulated"].astype(str).str.lower() == "true"


def _without_simulated(df: pd.DataFrame, include_simulated: bool) -> pd.DataFrame:
    if include_simulated or df.empty or "simulated" not in df.columns:
        return df
    return df[~_is_simulated(df)].reset_index(drop=True)


def load_registry(include_simulated: bool = False) -> pd.DataFrame:
    if not registry_exists():
        return pd.DataFrame()
    df = pd.read_csv(STUDY_REGISTRY_CSV)
    if include_simulated and DEMO_REGISTRY_CSV.exists() and DEMO_REGISTRY_CSV.stat().st_size > 0:
        df = pd.concat([df, pd.read_csv(DEMO_REGISTRY_CSV)], ignore_index=True)
    # The accepted name, so a synonym (Crassostrea gigas) filters with its species.
    df.insert(1, "species_accepted", df["species"].map(canonical_species))
    return df


def canonical_species(name) -> str:
    resolved = resolve_species(name) if isinstance(name, str) else None
    return resolved.scientific_name if resolved else (name if isinstance(name, str) else "")


def species_label(taxid) -> str:
    term = species_term_for_taxid(taxid)
    return term["scientific_name"] if term else f"taxid {taxid}"


def load_species_summary(include_simulated: bool = False) -> pd.DataFrame:
    """Per-species studies, crosswalk and evidence counts (reporting.species)."""
    studies = [load_yaml(p) for p in study_files(include_demo=include_simulated)]
    rows = species_summary(studies, load_evidence(include_simulated) if evidence_exists() else None)
    cols = ["species", "label", "crosswalk_route", "crosswalk_genes", "n_studies_real",
            "n_studies_real_harmonized", "n_evidence_records_real", "n_studies_simulated",
            "n_evidence_records_simulated", "complex"]
    return pd.DataFrame(rows)[cols]


def load_evidence(include_simulated: bool = False) -> pd.DataFrame:
    if not evidence_exists():
        return pd.DataFrame()
    return _without_simulated(load_evidence_table(), include_simulated)


def load_candidates(phenotype: str | None = None, include_simulated: bool = False) -> pd.DataFrame:
    """Candidates computed live from the evidence table.

    Real and simulated evidence are never pooled together, so dropping the
    simulated candidates afterwards leaves the real ranking unchanged.
    """
    evidence_df = load_evidence(include_simulated=True)
    if len(evidence_df) == 0:
        return pd.DataFrame()
    frames = []
    feature_types = sorted(evidence_df["feature_type"].dropna().unique())
    for ftype in feature_types:
        meta_df = run_meta_analysis(phenotype=phenotype, feature_type=ftype)
        if len(meta_df):
            frames.append(build_candidates(meta_df, evidence_df))
    if not frames:
        return pd.DataFrame()
    candidates = pd.concat(frames, ignore_index=True).sort_values(["tier", "score"], ascending=[True, False])
    return _without_simulated(candidates, include_simulated)


def phenotype_labels() -> dict:
    return {tid: term["label"] for tid, term in load_vocab("phenotype_ontology").items()}
