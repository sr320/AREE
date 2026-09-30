"""The dashboard builder must summarize whatever the pipeline has produced, keep
real and simulated evidence apart, and say plainly when outputs are missing."""
import json

import pytest
from click.testing import CliRunner
from conftest import ALL_DEMO_STUDY_IDS

import reporting.dashboard as dashboard
from aree.cli import main
from reporting.dashboard import build_dashboard_data

REQUIRED_TOP_LEVEL = {
    "generated_at", "repo_url", "commit", "artifacts_present", "kpis", "studies",
    "comparisons", "evidence", "candidates", "pools", "activity", "tests",
    "implementation_status", "vocab",
}


def _point_at_reports(monkeypatch, reports_dir):
    monkeypatch.setattr(dashboard, "EVIDENCE_TABLE", reports_dir / "evidence" / "evidence_table.tsv")
    monkeypatch.setattr(dashboard, "CANDIDATES_TABLE", reports_dir / "evidence_cards" / "candidates.tsv")
    monkeypatch.setattr(dashboard, "META_ANALYSIS_DIR", reports_dir / "meta_analysis")
    monkeypatch.setattr(dashboard, "MANIFESTS_DIR", reports_dir / "manifests")
    monkeypatch.setattr(dashboard, "TOP_CANDIDATES_MD", reports_dir / "top_candidates_summary.md")


def test_builder_reports_missing_outputs_instead_of_failing(tmp_path, monkeypatch):
    _point_at_reports(monkeypatch, tmp_path / "empty_reports")
    out = tmp_path / "site" / "dashboard" / "data.json"

    data = build_dashboard_data(out)

    assert out.exists()
    assert set(data) >= REQUIRED_TOP_LEVEL
    present = data["artifacts_present"]
    assert present["evidence_table"] is False
    assert present["candidates"] is False
    # Registry-derived counts still come through; pipeline-derived counts are zero.
    assert data["kpis"]["n_studies_simulated"] >= 6
    assert data["kpis"]["n_studies_real"] >= 1
    assert data["kpis"]["n_evidence_records_real"] == 0
    assert data["candidates"]["top_real"] == []
    assert all(not c["harmonized"] for c in data["comparisons"])
    # Nothing fabricated: no cards were copied when there were none to copy.
    assert not (out.parent / "cards").exists()


def test_builder_summarizes_a_demo_pipeline_run(isolated_reports, isolated_registry, tmp_path, monkeypatch):
    runner = CliRunner()
    for study_id in ALL_DEMO_STUDY_IDS:
        assert runner.invoke(main, ["register-study", f"registry/studies/{study_id}.yaml"]).exit_code == 0
        assert runner.invoke(main, ["harmonize", "--study", study_id, "--date", "2026-01-01"]).exit_code == 0
    assert runner.invoke(main, ["meta-analyze", "--feature-type", "gene"]).exit_code == 0
    assert runner.invoke(main, ["build-evidence-cards"]).exit_code == 0

    reports_dir = isolated_reports["evidence_table_path"].parent.parent
    _point_at_reports(monkeypatch, reports_dir)
    out = tmp_path / "site" / "dashboard" / "data.json"

    data = build_dashboard_data(out)

    assert data["artifacts_present"]["evidence_table"] is True
    assert data["artifacts_present"]["candidates"] is True
    kpis = data["kpis"]
    assert kpis["n_evidence_records_simulated"] > 0
    # The demo run harmonized no real study, so real evidence must read as zero,
    # not be inflated by the simulated records.
    assert kpis["n_evidence_records_real"] == 0
    assert kpis["n_candidates_real_high_priority"] == 0
    assert all(r["simulated"] for r in data["candidates"]["top_simulated"])
    assert data["candidates"]["top_real"] == []

    by_study = {r["study_id"]: r for r in data["evidence"]["by_study"]}
    assert set(ALL_DEMO_STUDY_IDS) <= set(by_study)
    assert all(by_study[s]["simulated"] for s in ALL_DEMO_STUDY_IDS)

    # Every demo comparison is marked harmonized and its record count matches.
    demo_comparisons = [c for c in data["comparisons"] if c["simulated"]]
    assert demo_comparisons and all(c["harmonized"] for c in demo_comparisons)
    assert sum(c["n_evidence_records"] for c in demo_comparisons) == kpis["n_evidence_records_simulated"]

    # The Pages workflow only runs `meta-analyze --feature-type gene`, which
    # writes the all-phenotypes file, so pools must be read from it too.
    meta_files = sorted(p.name for p in (reports_dir / "meta_analysis").glob("*_meta_analysis.tsv"))
    assert meta_files == ["all_phenotypes_gene_meta_analysis.tsv"]
    assert data["pools"], "no pools summarized from the all-phenotypes meta-analysis file"
    assert all(p["simulated"] for p in data["pools"])
    assert {p["phenotype"] for p in data["pools"]} >= {"thermal_tolerance", "larval_viability"}
    assert kpis["n_real_pools_k2plus"] == 0

    # Cards for the top simulated candidates were copied beside the JSON with
    # Quarto front matter, and the JSON links to them by site-relative path.
    linked = [r["card_page"] for r in data["candidates"]["top_simulated"] if r["card_page"]]
    assert linked
    for page in linked:
        card = out.parent / page.replace(".html", ".md")
        assert card.exists()
        text = card.read_text()
        assert text.startswith("---\ntitle:")
        assert "\nsidebar: false\n" in text.split("\n---\n", 1)[0]

    # The JSON round-trips (no NaN, no numpy scalars).
    reloaded = json.loads(out.read_text())
    assert reloaded["kpis"] == kpis


def test_cli_build_dashboard_writes_json(tmp_path, monkeypatch):
    _point_at_reports(monkeypatch, tmp_path / "empty_reports")
    out = tmp_path / "data.json"
    result = CliRunner().invoke(main, ["build-dashboard", "--out", str(out)])
    assert result.exit_code == 0, result.output
    assert "Pipeline outputs not found" in result.output
    assert json.loads(out.read_text())["kpis"]["n_studies_simulated"] >= 6


@pytest.mark.parametrize("citation,expected", [
    ("Calla B, Thompson NF, Burge CA (2026). Population-specific ...", "Calla B, Thompson NF, Burge CA (2026)"),
    ("", ""),
    (None, ""),
])
def test_short_citation(citation, expected):
    assert dashboard._short_citation(citation) == expected
