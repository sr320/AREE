"""End-to-end smoke test: the full demo pipeline (validate -> register -> harmonize ->
meta-analyze -> build-evidence-cards) must run successfully on every simulated demo study."""
from click.testing import CliRunner
from conftest import ALL_DEMO_STUDY_IDS

from aree.cli import main


def test_full_demo_pipeline_via_cli(isolated_reports, isolated_registry, tmp_path, monkeypatch):
    runner = CliRunner()

    for study_id in ALL_DEMO_STUDY_IDS:
        result = runner.invoke(main, ["validate-study", f"registry/demo_studies/{study_id}.yaml"])
        assert result.exit_code == 0, result.output

    for study_id in ALL_DEMO_STUDY_IDS:
        result = runner.invoke(main, ["register-study", f"registry/demo_studies/{study_id}.yaml"])
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
    runner.invoke(main, ["register-study", "registry/demo_studies/GIGAS_HEAT01.yaml"])
    result = runner.invoke(main, ["register-study", "registry/demo_studies/GIGAS_HEAT01.yaml"])
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


def test_metabolomics_handoff_carries_standard_error_and_method(tmp_path):
    """The metabolomics raw mode writes lfcSE and a per-row analysis_method;
    both must reach the evidence record. A processed table without them (the
    demo's) keeps the old defaults."""
    import pandas as pd

    from harmonize.core import harmonize_processed_table

    out = tmp_path / "GIGAS_GROW06_hypoxia_vs_control_metabolite_features_standardized.tsv"
    pd.DataFrame({
        "feature_id": ["GM_alanine", "GM_glycine"],
        "putative_metabolite_name": ["alanine", "glycine"],
        "annotation_confidence_level": [1, 1],
        "metabolite_class": ["amino_acid", "amino_acid"],
        "log2FC": [1.98, -1.61],
        "pvalue": [3e-10, 5e-9],
        "padj": [3e-9, 2e-8],
        "lfcSE": [0.032, 0.041],
        "analysis_method": ["targeted_MS_limma_moderated_t"] * 2,
    }).to_csv(out, sep="\t", index=False)
    df = harmonize_processed_table("GIGAS_GROW06", out, date_generated="2026-01-01",
                                   comparison_id="hypoxia_vs_control")
    assert df["standard_error"].tolist() == [0.032, 0.041]
    assert set(df["analysis_method"]) == {"targeted_MS_limma_moderated_t"}

    legacy = out.with_name("legacy.tsv")
    pd.read_csv(out, sep="\t").drop(columns=["lfcSE", "analysis_method"]).to_csv(legacy, sep="\t", index=False)
    df = harmonize_processed_table("GIGAS_GROW06", legacy, date_generated="2026-01-01",
                                   comparison_id="hypoxia_vs_control")
    assert df["standard_error"].isna().all()
    assert set(df["analysis_method"]) == {"untargeted_LCMS"}


def test_proteomics_handoff_carries_standard_error_method_and_legacy_ids(isolated_reports, tmp_path):
    """The proteomics raw mode writes lfcSE and analysis_method; both must reach
    the evidence record. oyster_v9 CGI_ models take the legacy-locus-tag route
    (unresolved against the demo crosswalk, which has no such column), not the
    UniProt one. The demo's processed table keeps the old defaults."""
    import pandas as pd

    from harmonize.core import harmonize_processed_table
    from harmonize.proteomics import protein_id_type

    assert protein_id_type("CGI_10001851") == "legacy_locus_tag"
    assert protein_id_type("CGI_10000041,CGI_10008574") == "legacy_locus_tag"
    assert protein_id_type("K1QN12") == "uniprot_accession"

    out = tmp_path / "GIGAS_SAL04_low_salinity_vs_control_protein_abundance_standardized.tsv"
    pd.DataFrame({
        "protein_accession": ["K1QN12", "CGI_10001851"],
        "gene_symbol": ["hsp70", "CGI_10001851"],
        "log2FC": [1.4, -0.8],
        "lfcSE": [0.21, 0.33],
        "pvalue": [0.001, 0.02],
        "padj": [0.01, 0.08],
        "missingness_percent": [0.0, 12.5],
        "analysis_method": ["spectral_count_limma_moderated_t"] * 2,
    }).to_csv(out, sep="\t", index=False)
    df = harmonize_processed_table("GIGAS_SAL04", out, date_generated="2026-01-01",
                                   comparison_id="low_salinity_vs_control")
    assert df["standard_error"].tolist() == [0.21, 0.33]
    assert set(df["analysis_method"]) == {"spectral_count_limma_moderated_t"}
    assert df.set_index("feature_id_original").loc["CGI_10001851", "mapping_confidence"] == "unresolved"
    assert df.set_index("feature_id_original").loc["K1QN12", "mapping_confidence"] != "unresolved"

    legacy = out.with_name("legacy.tsv")
    pd.read_csv(out, sep="\t").drop(columns=["lfcSE", "analysis_method"]).to_csv(legacy, sep="\t", index=False)
    df = harmonize_processed_table("GIGAS_SAL04", legacy, date_generated="2026-01-01",
                                   comparison_id="low_salinity_vs_control")
    assert df["standard_error"].isna().all()
    assert set(df["analysis_method"]) == {"MSstats"}


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
