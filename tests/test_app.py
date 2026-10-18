"""The Streamlit app runs headless without errors and shows every registered species. WP8.

The app reads the evidence table through `data_access` and `harmonize.core`;
both are pointed at a temporary demo run so the test never depends on what
happens to be in reports/.
"""
from __future__ import annotations

import sys

import pytest

pytest.importorskip("streamlit")
from conftest import harmonize_all_demo_studies  # noqa: E402
from streamlit.testing.v1 import AppTest  # noqa: E402

from common import REPO_ROOT, load_vocab  # noqa: E402

APP = REPO_ROOT / "app" / "main.py"


@pytest.fixture
def app_on(isolated_reports, monkeypatch):
    sys.path.insert(0, str(APP.parent))
    import data_access

    monkeypatch.setattr(data_access, "EVIDENCE_TABLE_PATH", isolated_reports["evidence_table_path"])
    yield isolated_reports
    sys.path.remove(str(APP.parent))


def _run():
    at = AppTest.from_file(str(APP), default_timeout=180).run()
    assert not at.exception, at.exception
    return at


def test_app_without_evidence_still_lists_every_species(app_on):
    at = _run()
    species = at.dataframe[0].value
    assert set(species["species"]) == {t["scientific_name"] for t in load_vocab("species").values()}
    assert {"ncbi_gene", "gff_homology"} <= set(species["crosswalk_route"])
    assert species["n_evidence_records_real"].sum() == 0


def test_app_hides_simulated_demo_data_until_asked(app_on):
    harmonize_all_demo_studies()
    at = _run()
    table = at.dataframe[0].value.set_index("species")
    assert table["n_evidence_records_simulated"].sum() == 0
    assert not any("feature_id_standardized" in d.value.columns for d in at.dataframe)


def test_app_on_the_demo_filters_by_species(app_on):
    harmonize_all_demo_studies()
    at = _run()
    at.sidebar.checkbox(key="include_simulated").check().run()
    assert not at.exception, at.exception
    table = at.dataframe[0].value.set_index("species")
    assert table.loc["Mytilus galloprovincialis", "n_evidence_records_simulated"] > 0
    assert table.loc["Mytilus galloprovincialis", "complex"] == "mytilus_edulis_complex"
    # Filter the evidence tab to the demo mussel: every row shown is that species.
    # Selected by key: the studies tab has no filter when the registry index is
    # absent (as in CI, where tests run before it is built), so position lies.
    at.multiselect(key="evidence_species").set_value(["Mytilus galloprovincialis"]).run()
    shown = [d.value for d in at.dataframe
             if {"species_accepted", "feature_id_standardized"} <= set(d.value.columns)]
    assert shown and set(shown[0]["species_accepted"]) == {"Mytilus galloprovincialis"}


def test_candidate_loader_orders_stronger_tiers_before_emerging(app_on, monkeypatch):
    import data_access
    import pandas as pd

    monkeypatch.setattr(data_access, "load_evidence", lambda **kwargs: pd.DataFrame({
        "feature_type": ["gene"],
    }))
    monkeypatch.setattr(data_access, "run_meta_analysis", lambda **kwargs: pd.DataFrame({"x": [1]}))
    monkeypatch.setattr(data_access, "build_candidates", lambda *args: pd.DataFrame({
        "tier": ["emerging", "multi_omics_convergence", "high_priority_cross_study"],
        "score": [90.0, 70.0, 60.0], "simulated": [False, False, False],
    }))
    assert data_access.load_candidates()["tier"].tolist() == [
        "high_priority_cross_study", "multi_omics_convergence", "emerging",
    ]
