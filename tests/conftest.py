import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from common import EVIDENCE_TABLE_PATH, study_ids  # noqa: E402


def _stat(path: Path):
    return path.stat().st_mtime_ns if path.exists() else None


@pytest.fixture(autouse=True)
def _real_evidence_table_is_untouched():
    """Fail any test that writes the real reports/evidence/evidence_table.tsv.

    A harmonize call without `isolated_reports` upserts its fixture rows into the
    working evidence table, where they silently replace real or demo records
    until someone re-harmonizes (a stray GIGAS_GROW06 fixture once broke the
    full meta-analysis this way).
    """
    before = _stat(EVIDENCE_TABLE_PATH)
    yield
    assert _stat(EVIDENCE_TABLE_PATH) == before, (
        "this test modified the real evidence table; request the `isolated_reports` fixture"
    )


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


# Every simulated study in the registry, oyster and mussel alike. Derived, not
# listed, so a new demo study is exercised by every demo-pipeline test.
ALL_DEMO_STUDY_IDS = study_ids(simulated=True)


# The demo now holds two species that share feature identifiers on purpose
# (MYTGAL_HEAT01 reuses the oyster demo's IDs). Tests about oyster results must
# select the oyster partition explicitly rather than take the first row for a gene.
OYSTER_TAXID = 29159


def oyster(df):
    """Rows of a meta-analysis, candidate or evidence table in the oyster partition."""
    return df[df["species_taxid"].astype(int) == OYSTER_TAXID]


def oyster_card(index, feature_id):
    """The oyster evidence-card index entry for a feature."""
    return next(
        e for e in index
        if e["feature_id_standardized"] == feature_id and int(e["species_taxid"]) == OYSTER_TAXID
    )


def harmonize_all_demo_studies(date_generated="2026-01-01"):
    from harmonize.core import harmonize_study
    for study_id in ALL_DEMO_STUDY_IDS:
        harmonize_study(study_id, date_generated=date_generated)
