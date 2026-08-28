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


# --------------------------------------------------------------------------- #
# Workflow output -> harmonize handoff
# --------------------------------------------------------------------------- #


def test_harmonize_accepts_an_explicit_comparison(tmp_path):
    """The Nextflow workflow names its output for the pipeline stage that made
    it, not for any path in the registry, so the comparison cannot always be
    inferred from the filename. --comparison is how workflow output connects."""
    import pandas as pd

    from harmonize.core import harmonize_processed_table

    # Same columns STANDARDIZE_OUTPUT emits, under a workflow-style filename
    # that matches no declared results_file.
    out = tmp_path / "GIGAS_HEAT01_acute_heat_vs_control_dge_standardized.tsv"
    pd.DataFrame({
        "gene_id": ["LOC105333935", "LOC105333100"],
        "baseMean": [4210.5, 3105.2],
        "log2FoldChange": [2.4, 1.8],
        "lfcSE": [0.31, 0.29],
        "stat": [7.74, 6.21],
        "pvalue": [9.8e-15, 5.3e-10],
        "padj": [0.0008, 0.003],
    }).to_csv(out, sep="\t", index=False)

    df = harmonize_processed_table(
        "GIGAS_HEAT01", out, date_generated="2026-01-01",
        comparison_id="acute_heat_vs_control",
    )
    assert len(df) == 2
    assert set(df["comparison_id"]) == {"acute_heat_vs_control"}
    # lfcSE survives the handoff — this is what makes a raw-reanalysis study poolable.
    assert df["standard_error"].notna().all()


def test_harmonize_rejects_an_unknown_comparison(tmp_path):
    import pandas as pd
    import pytest

    from harmonize.core import harmonize_processed_table

    out = tmp_path / "whatever.tsv"
    pd.DataFrame({"gene_id": ["LOC105333935"], "log2FoldChange": [1.0],
                  "lfcSE": [0.2], "pvalue": [0.01], "padj": [0.05],
                  "baseMean": [100.0], "stat": [5.0]}).to_csv(out, sep="\t", index=False)

    with pytest.raises(ValueError, match="has no comparison"):
        harmonize_processed_table("GIGAS_HEAT01", out, date_generated="2026-01-01",
                                  comparison_id="no_such_comparison")


def test_unmatched_filename_without_comparison_says_what_to_do(tmp_path):
    import pandas as pd
    import pytest

    from harmonize.core import harmonize_processed_table

    out = tmp_path / "mystery_output.tsv"
    pd.DataFrame({"gene_id": ["LOC105333935"], "log2FoldChange": [1.0],
                  "lfcSE": [0.2], "pvalue": [0.01], "padj": [0.05],
                  "baseMean": [100.0], "stat": [5.0]}).to_csv(out, sep="\t", index=False)

    with pytest.raises(ValueError, match="--comparison"):
        harmonize_processed_table("GIGAS_HEAT01", out, date_generated="2026-01-01")
