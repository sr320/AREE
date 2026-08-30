import pandas as pd
from conftest import harmonize_all_demo_studies

from harmonize.core import load_evidence_table
from meta_analysis.run import run_meta_analysis
from prioritize.rank import build_candidates
from prioritize.scoring import candidate_score, compute_components


def test_candidate_score_is_pure_and_reproducible():
    components = {
        "n_studies_score": 0.4, "sample_size_score": 0.5, "effect_magnitude_score": 0.6,
        "significance_score": 0.7, "direction_consistency_score": 1.0,
        "phenotype_relevance_score": 1.0, "context_breadth_score": 0.3,
        "assay_diversity_score": 0.3, "mapping_confidence_score": 1.0,
        "quality_score": 0.8, "heterogeneity_penalty": 0.1,
    }
    score1 = candidate_score(components)
    score2 = candidate_score(components)
    assert score1 == score2


def test_candidate_score_penalizes_heterogeneity():
    base = {
        "n_studies_score": 1.0, "sample_size_score": 1.0, "effect_magnitude_score": 1.0,
        "significance_score": 1.0, "direction_consistency_score": 1.0,
        "phenotype_relevance_score": 1.0, "context_breadth_score": 1.0,
        "assay_diversity_score": 1.0, "mapping_confidence_score": 1.0, "quality_score": 1.0,
    }
    low_het = candidate_score({**base, "heterogeneity_penalty": 0.0})
    high_het = candidate_score({**base, "heterogeneity_penalty": 1.0})
    assert high_het < low_het


def test_compute_components_bounds_are_0_to_1():
    meta_row = {
        "k_studies": 10, "total_sample_size": 10000, "pooled_effect": 50.0,
        "p_value": 1e-300, "direction_consistency": 1.0, "phenotype": "survival",
        "distinct_tissues": 10, "distinct_life_stages": 10, "i_squared": 100.0,
        "n_distinct_assays": 10,
    }
    components = compute_components(meta_row, ["exact"], [])
    for key, value in components.items():
        assert 0.0 <= value <= 1.0, f"{key}={value} out of bounds"


def test_high_priority_gate_blocks_low_direction_consistency_even_with_two_studies(isolated_reports):
    harmonize_all_demo_studies()
    meta_df = run_meta_analysis(phenotype="larval_viability", feature_type="gene")
    evidence_df = load_evidence_table()
    candidates = build_candidates(meta_df, evidence_df)

    sod1 = candidates[candidates["feature_id_standardized"] == "LOC105331241"].iloc[0]
    assert sod1["k_studies"] == 2
    assert sod1["tier"] != "high_priority_cross_study", (
        "sod1 has 2 studies but conflicting direction and must not be promoted to high_priority"
    )


def test_hsp70_reaches_high_priority_tier(isolated_reports):
    harmonize_all_demo_studies()
    meta_df = run_meta_analysis(phenotype="larval_viability", feature_type="gene")
    evidence_df = load_evidence_table()
    candidates = build_candidates(meta_df, evidence_df)

    hsp70 = candidates[candidates["feature_id_standardized"] == "LOC105333935"].iloc[0]
    assert hsp70["tier"] == "high_priority_cross_study"


def test_assay_diversity_never_crosses_origin_or_species_partitions():
    meta_df = pd.DataFrame([{
        "feature_id_standardized": "GENE1",
        "phenotype": "survival",
        "feature_type": "gene",
        "simulated": False,
        "species_taxid": 29159,
        "k_studies": 1,
        "studies": "REAL_GENE",
        "n_evidence_records": 1,
        "total_sample_size": 12,
        "pooled_effect": 1.0,
        "pooled_se": 0.2,
        "ci_lower": 0.6,
        "ci_upper": 1.4,
        "z": 5.0,
        "p_value": 1e-4,
        "q_statistic": 0.0,
        "i_squared": 0.0,
        "tau_squared": 0.0,
        "direction_consistency": 1.0,
        "distinct_tissues": 1,
        "distinct_life_stages": 1,
        "distinct_stressors": "temperature",
        "contributing_evidence_ids": "real-gene",
    }])
    evidence_df = pd.DataFrame([
        {
            "evidence_id": "real-gene", "feature_id_standardized": "GENE1",
            "phenotype": "survival", "feature_type": "gene", "simulated": False,
            "species_taxid": 29159, "mapping_confidence": "exact", "quality_flags": [],
        },
        {
            "evidence_id": "sim-protein", "feature_id_standardized": "GENE1",
            "phenotype": "survival", "feature_type": "protein", "simulated": True,
            "species_taxid": 29159, "mapping_confidence": "inferred",
            "quality_flags": ["processed_only"],
        },
        {
            "evidence_id": "other-species-protein", "feature_id_standardized": "GENE1",
            "phenotype": "survival", "feature_type": "protein", "simulated": False,
            "species_taxid": 99999, "mapping_confidence": "exact", "quality_flags": [],
        },
    ])

    candidate = build_candidates(meta_df, evidence_df).iloc[0]
    assert candidate["n_distinct_assays"] == 1
    assert candidate["mapping_confidences"] == "exact"
    assert candidate["quality_flags_union"] == ""
