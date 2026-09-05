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


# --------------------------------------------------------------------------- #
# Numerical precision and multiple testing
# --------------------------------------------------------------------------- #


def test_pooled_p_value_does_not_underflow_to_zero():
    """z = 20 gives p ~ 5e-89; `2 * (1 - cdf)` returned exactly 0.0 here and
    printed as 'p = 0' on real evidence cards."""
    result = dersimonian_laird([4.0], [0.2])
    assert result.z == 20.0
    assert 0.0 < result.p_value < 1e-80


def test_benjamini_hochberg_matches_reference_values():
    from meta_analysis.pooling import benjamini_hochberg

    adjusted = benjamini_hochberg([0.01, 0.04, 0.03, 0.20])
    # ranks: 0.01 -> 0.01*4/1 = 0.04; 0.03 -> 0.03*4/2 = 0.06; 0.04 -> 0.04*4/3 = 0.0533,
    # then the step-up makes 0.03's value min(0.06, 0.0533) = 0.0533; 0.20 -> 0.20.
    assert [round(a, 4) for a in adjusted] == [0.04, 0.0533, 0.0533, 0.2]
    assert (adjusted <= 1.0).all()


def test_benjamini_hochberg_leaves_nan_out_of_the_family():
    import math

    from meta_analysis.pooling import benjamini_hochberg

    adjusted = benjamini_hochberg([0.02, float("nan"), 0.5])
    assert math.isnan(adjusted[1])
    # n = 2, not 3: 0.02 * 2 / 1 = 0.04
    assert math.isclose(adjusted[0], 0.04)


def test_adjusted_p_is_computed_within_family_and_invariant_to_filters(isolated_reports):
    """The family is (phenotype, feature_type, simulated, species_taxid), which is
    exactly what the CLI can filter on, so the same feature gets the same adjusted
    p whether the run covered one phenotype or every phenotype."""
    harmonize_all_demo_studies()
    everything = run_meta_analysis(feature_type="gene")
    one_phenotype = run_meta_analysis(phenotype="larval_viability", feature_type="gene")

    assert "adjusted_p_value" in everything.columns
    assert (everything["adjusted_p_value"] >= everything["p_value"] - 1e-12).all()

    hsp70_all = everything[
        (everything["feature_id_standardized"] == "LOC105333935")
        & (everything["phenotype"] == "larval_viability")
    ].iloc[0]
    hsp70_one = one_phenotype[one_phenotype["feature_id_standardized"] == "LOC105333935"].iloc[0]
    assert hsp70_all["adjusted_p_value"] == hsp70_one["adjusted_p_value"]
    assert hsp70_all["n_tests_in_family"] == hsp70_one["n_tests_in_family"] == len(one_phenotype)


def test_numeric_gene_ids_survive_as_text(isolated_reports):
    """Real NCBI GeneIDs are all digits. If the meta-analysis reads them as
    integers while the evidence loader keeps text, no candidate can find its
    own evidence records."""
    harmonize_all_demo_studies()
    evidence = load_evidence_table()
    row = evidence.iloc[0].copy()
    row["feature_id_standardized"] = "105317636"
    row["feature_id_original"] = "105317636"
    row["evidence_id"] = "numeric-id-row"
    row["standard_error"] = 0.2
    row["p_value"] = 0.01
    pd.DataFrame([row]).to_csv(isolated_reports["evidence_table_path"], sep="\t", index=False)

    result = run_meta_analysis(feature_type=row["feature_type"])
    assert result["feature_id_standardized"].tolist() == ["105317636"]
    # pandas 2 keeps text as object dtype, pandas 3 as StringDtype; the contract
    # is that the value is text, not which dtype holds it.
    assert isinstance(result["feature_id_standardized"].iloc[0], str)


def test_vectorized_standard_errors_match_the_scalar_definition():
    """`effective_standard_errors` is the whole-table twin of
    `effective_standard_error`; they must agree on every edge case."""
    import numpy as np

    from meta_analysis.effect_sizes import effective_standard_error
    from meta_analysis.run import effective_standard_errors

    rows = pd.DataFrame({
        "effect_size":    [1.0, 1.0,  -2.0, 0.0,  1.5,  1.5, None, 0.8,  3.0],
        "standard_error": [0.2, None, None, None, 0.0,  -1., 0.3,  None, None],
        "p_value":        [0.5, 0.05, 1e-9, 0.01, 0.02, 0.5, 0.01, None, 0.0],
    })
    vectorized = effective_standard_errors(rows)
    for i, row in rows.iterrows():
        scalar = effective_standard_error(row["effect_size"], row["standard_error"], row["p_value"])
        # The scalar helper returns None for missing inputs but propagates NaN
        # when the p-value itself is NaN; both mean "no usable SE".
        if scalar is None or scalar != scalar:
            assert np.isnan(vectorized[i]), f"row {i}: scalar missing, vectorized {vectorized[i]}"
        else:
            assert np.isclose(vectorized[i], scalar), f"row {i}: {vectorized[i]} vs {scalar}"
