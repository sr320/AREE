"""End-to-end smoke test: the full demo pipeline (validate -> register -> harmonize ->
meta-analyze -> build-evidence-cards) must run successfully on all 6 demo studies."""
from click.testing import CliRunner
from conftest import ALL_DEMO_STUDY_IDS

from aree.cli import main


def test_full_demo_pipeline_via_cli(isolated_reports, isolated_registry, tmp_path, monkeypatch):
    runner = CliRunner()

    for study_id in ALL_DEMO_STUDY_IDS:
        result = runner.invoke(main, ["validate-study", f"registry/studies/{study_id}.yaml"])
        assert result.exit_code == 0, result.output

    for study_id in ALL_DEMO_STUDY_IDS:
        result = runner.invoke(main, ["register-study", f"registry/studies/{study_id}.yaml"])
        assert result.exit_code == 0, result.output

    for study_id in ALL_DEMO_STUDY_IDS:
        result = runner.invoke(main, ["harmonize", "--study", study_id, "--date", "2026-01-01"])
        assert result.exit_code == 0, result.output

    result = runner.invoke(main, ["meta-analyze", "--feature-type", "gene"])
    assert result.exit_code == 0, result.output

    result = runner.invoke(main, ["build-evidence-cards"])
    assert result.exit_code == 0, result.output
    assert "Wrote" in result.output


def test_duplicate_registration_exits_nonzero_without_update_flag(isolated_registry):
    runner = CliRunner()
    runner.invoke(main, ["register-study", "registry/studies/GIGAS_HEAT01.yaml"])
    result = runner.invoke(main, ["register-study", "registry/studies/GIGAS_HEAT01.yaml"])
    assert result.exit_code != 0
