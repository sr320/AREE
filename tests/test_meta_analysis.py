import math

from conftest import harmonize_all_demo_studies
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
    assert "unresolved" not in set(
        conf for confs in result["mapping_confidences"] for conf in confs.split("|")
    )
