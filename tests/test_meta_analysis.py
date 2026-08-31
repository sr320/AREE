import math

import pandas as pd
import pytest
from conftest import harmonize_all_demo_studies

from harmonize.core import load_evidence_table
from meta_analysis.effect_sizes import approximate_se
from meta_analysis.pooling import dersimonian_laird
from meta_analysis.run import run_meta_analysis


def test_approximate_se_recovers_z_relationship():
    # p = 0.05 two-sided corresponds to z ~= 1.96
    se = approximate_se(effect_size=1.0, p_value=0.05)
    assert math.isclose(se, 1.0 / 1.959963985, rel_tol=1e-3)


def test_dersimonian_laird_single_study_is_fixed_effect():
    result = dersimonian_laird([1.5], [0.3])
    assert result.k == 1
    assert math.isclose(result.pooled_effect, 1.5)
    assert result.i_squared == 0.0
    assert result.tau_squared == 0.0


def test_dersimonian_laird_two_consistent_studies_pool_toward_shared_effect():
    result = dersimonian_laird([1.0, 1.2], [0.2, 0.25])
    assert result.k == 2
    assert 1.0 <= result.pooled_effect <= 1.2
    assert result.p_value < 0.05


def test_dersimonian_laird_conflicting_studies_show_high_heterogeneity():
    result = dersimonian_laird([1.0, -1.0], [0.2, 0.2])
    assert result.i_squared > 50.0


def test_meta_analysis_on_demo_data_pools_hsp70_across_studies(isolated_reports):
    harmonize_all_demo_studies()
    result = run_meta_analysis(phenotype="larval_viability", feature_type="gene")
    hsp70 = result[result["feature_id_standardized"] == "LOC105333935"].iloc[0]
    assert hsp70["k_studies"] == 2
    assert hsp70["direction_consistency"] == 1.0


def test_meta_analysis_flags_conflicting_direction_for_sod1(isolated_reports):
    harmonize_all_demo_studies()
    result = run_meta_analysis(phenotype="larval_viability", feature_type="gene")
    sod1 = result[result["feature_id_standardized"] == "LOC105331241"].iloc[0]
    assert sod1["k_studies"] == 2
    assert sod1["direction_consistency"] < 0.7
    assert sod1["i_squared"] > 50.0


def test_meta_analysis_excludes_unresolved_identifiers(isolated_reports):
    harmonize_all_demo_studies()
    result = run_meta_analysis(feature_type="gene")
    assert "unresolved" not in {
        conf for confs in result["mapping_confidences"] for conf in confs.split("|")
    }


def test_unpoolable_study_does_not_count_as_replication(isolated_reports):
    harmonize_all_demo_studies()
    evidence = load_evidence_table()
    source = evidence[
        (evidence["feature_id_standardized"] == "LOC105331241")
        & (evidence["phenotype"] == "larval_viability")
    ].iloc[0]

    poolable = source.copy()
    poolable["evidence_id"] = "poolable"
    poolable["study_id"] = "POOLABLE_STUDY"
    poolable["comparison_id"] = "poolable_comparison"
    poolable["standard_error"] = 0.2
    poolable["p_value"] = 0.01

    unavailable = source.copy()
    unavailable["evidence_id"] = "unavailable"
    unavailable["study_id"] = "UNAVAILABLE_STUDY"
    unavailable["comparison_id"] = "unavailable_comparison"
    unavailable["standard_error"] = None
    unavailable["p_value"] = None

    path = isolated_reports["evidence_table_path"]
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([poolable, unavailable]).to_csv(path, sep="\t", index=False)

    result = run_meta_analysis(phenotype="larval_viability", feature_type="gene").iloc[0]
    assert result["k_studies"] == 1
    assert result["studies"] == "POOLABLE_STUDY"
    assert result["n_evidence_records"] == 1
    assert result["n_available_records"] == 2
    assert result["n_excluded_unpoolable"] == 1
    assert result["excluded_studies"] == "UNAVAILABLE_STUDY"


def test_multiple_comparisons_from_one_study_fail_closed(isolated_reports):
    harmonize_all_demo_studies()
    evidence = load_evidence_table()
    source = evidence[
        (evidence["feature_id_standardized"] == "LOC105331241")
        & (evidence["phenotype"] == "larval_viability")
    ].iloc[0]

    rows = []
    for comparison_id, effect in [("contrast_a", 1.0), ("contrast_b", 1.2)]:
        row = source.copy()
        row["evidence_id"] = comparison_id
        row["study_id"] = "SHARED_STUDY"
        row["comparison_id"] = comparison_id
        row["effect_size"] = effect
        row["standard_error"] = 0.2
        row["p_value"] = 0.01
        rows.append(row)

    path = isolated_reports["evidence_table_path"]
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)

    with pytest.raises(ValueError, match="multiple comparisons from one study"):
        run_meta_analysis(phenotype="larval_viability", feature_type="gene")


def test_lower_confidence_alias_is_not_pooled_as_a_second_effect(isolated_reports):
    harmonize_all_demo_studies()
    result = run_meta_analysis(phenotype="larval_viability", feature_type="gene")
    hsp70 = result[result["feature_id_standardized"] == "LOC105333935"].iloc[0]
    assert hsp70["n_excluded_duplicate_mappings"] == 1
    assert hsp70["n_excluded_unpoolable"] == 0
    assert hsp70["n_evidence_records"] == hsp70["k_studies"] == 2
