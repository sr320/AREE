import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))


@pytest.fixture
def isolated_reports(tmp_path, monkeypatch):
    """Redirect the evidence table / meta-analysis / evidence-card output paths
    to a temp directory so tests never mutate the real reports/ directory."""
    import harmonize.core as harmonize_core
    import meta_analysis.run as meta_run
    import reporting.evidence_cards as evidence_cards

    evidence_table_path = tmp_path / "evidence" / "evidence_table.tsv"
    meta_dir = tmp_path / "meta_analysis"
    cards_dir = tmp_path / "evidence_cards"
    manifests_dir = tmp_path / "manifests"

    monkeypatch.setattr(harmonize_core, "EVIDENCE_TABLE_PATH", evidence_table_path)
    monkeypatch.setattr(harmonize_core, "MANIFESTS_DIR", manifests_dir)
    monkeypatch.setattr(meta_run, "EVIDENCE_TABLE_PATH", evidence_table_path)
    monkeypatch.setattr(meta_run, "META_ANALYSIS_DIR", meta_dir)
    monkeypatch.setattr(evidence_cards, "EVIDENCE_CARDS_DIR", cards_dir)
    monkeypatch.setattr(evidence_cards, "REPORTS_DIR", tmp_path)

    return {
        "evidence_table_path": evidence_table_path, "meta_dir": meta_dir,
        "cards_dir": cards_dir, "manifests_dir": manifests_dir,
    }


@pytest.fixture
def isolated_registry(tmp_path, monkeypatch):
    """Redirect the study registry CSV to a temp file so tests never mutate
    the real registry/study_registry.csv."""
    import intake.registry as registry_mod

    registry_path = tmp_path / "study_registry.csv"
    monkeypatch.setattr(registry_mod, "STUDY_REGISTRY_CSV", registry_path)
    return registry_path


ALL_DEMO_STUDY_IDS = [
    "GIGAS_HEAT01", "GIGAS_OA02", "GIGAS_PATH03",
    "GIGAS_SAL04", "GIGAS_LARV05", "GIGAS_GROW06",
]


def harmonize_all_demo_studies(date_generated="2026-01-01"):
    from harmonize.core import harmonize_study
    for study_id in ALL_DEMO_STUDY_IDS:
        harmonize_study(study_id, date_generated=date_generated)
