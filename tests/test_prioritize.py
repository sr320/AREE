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
    assert candidate["n_supporting_layers"] == 1
    assert candidate["mapping_confidences"] == "exact"
    assert candidate["quality_flags_union"] == ""


def _two_study_meta_row(**overrides) -> dict:
    row = {
        "feature_id_standardized": "GENE_NULL", "phenotype": "disease_susceptibility",
        "feature_type": "gene", "simulated": False, "species_taxid": 29159,
        "k_studies": 2, "studies": "A|B", "n_evidence_records": 2, "total_sample_size": 17,
        "pooled_effect": 0.12, "pooled_se": 0.3, "ci_lower": -0.47, "ci_upper": 0.71, "z": 0.4,
        "p_value": 0.69, "adjusted_p_value": 0.93, "n_tests_in_family": 23094,
        "q_statistic": 0.1, "i_squared": 0.0, "tau_squared": 0.0,
        "direction_consistency": 1.0, "distinct_tissues": 1, "distinct_life_stages": 1,
        "distinct_stressors": "pathogen_challenge", "contributing_evidence_ids": "a|b",
    }
    row.update(overrides)
    return row


def _two_study_evidence() -> pd.DataFrame:
    return pd.DataFrame([
        {"evidence_id": "a", "study_id": "A", "feature_id_standardized": "GENE_NULL",
         "phenotype": "disease_susceptibility", "feature_type": "gene", "simulated": False,
         "species_taxid": 29159, "mapping_confidence": "exact", "quality_flags": []},
        {"evidence_id": "b", "study_id": "B", "feature_id_standardized": "GENE_NULL",
         "phenotype": "disease_susceptibility", "feature_type": "gene", "simulated": False,
         "species_taxid": 29159, "mapping_confidence": "exact", "quality_flags": []},
    ])


def test_high_priority_gate_requires_adjusted_significance():
    """Two genome-wide studies agree in sign for about half of all null genes.
    Replication + direction agreement without a significant pooled effect must
    not reach the top tier."""
    candidate = build_candidates(pd.DataFrame([_two_study_meta_row()]), _two_study_evidence()).iloc[0]
    assert candidate["k_studies"] == 2
    assert candidate["direction_consistency"] == 1.0
    assert candidate["tier"] == "emerging"
    assert not candidate["is_high_priority"]


def test_high_priority_gate_uses_adjusted_not_raw_p():
    """Nominally significant but not after BH: still not high priority."""
    row = _two_study_meta_row(p_value=0.004, adjusted_p_value=0.31)
    candidate = build_candidates(pd.DataFrame([row]), _two_study_evidence()).iloc[0]
    assert candidate["tier"] == "emerging"

    row = _two_study_meta_row(p_value=1e-6, adjusted_p_value=0.003, pooled_effect=1.4)
    candidate = build_candidates(pd.DataFrame([row]), _two_study_evidence()).iloc[0]
    assert candidate["tier"] == "high_priority_cross_study"


def test_candidates_are_ordered_strongest_tier_first(isolated_reports):
    from prioritize.rank import TIER_ORDER

    harmonize_all_demo_studies()
    meta_df = run_meta_analysis(feature_type="gene")
    candidates = build_candidates(meta_df, load_evidence_table())
    ranks = candidates["tier"].map({t: i for i, t in enumerate(TIER_ORDER)}).tolist()
    assert ranks == sorted(ranks)
    assert candidates["tier"].iloc[0] == "high_priority_cross_study"


def test_ranking_matches_candidates_to_numeric_gene_ids():
    """Candidate rows and evidence rows must agree on identifier type, or the
    forest data and mapping confidences silently come out empty."""
    meta = pd.DataFrame([_two_study_meta_row(feature_id_standardized="105317636")])
    evidence = _two_study_evidence().assign(feature_id_standardized="105317636")
    candidate = build_candidates(meta, evidence).iloc[0]
    assert candidate["mapping_confidences"] == "exact"


# --------------------------------------------------------------------------- #
# Multi-omics convergence: same phenotype, significant in each layer
# --------------------------------------------------------------------------- #


def _gene_candidate(**overrides) -> pd.DataFrame:
    row = _two_study_meta_row(
        feature_id_standardized="GENE_X", phenotype="thermal_tolerance", k_studies=1,
        studies="RNA", n_evidence_records=1, contributing_evidence_ids="rna",
        p_value=1e-6, adjusted_p_value=1e-4, pooled_effect=1.8,
    )
    row.update(overrides)
    return pd.DataFrame([row])


def _record(evidence_id, study_id, feature_type, phenotype, adjusted_p, feature="GENE_X") -> dict:
    return {
        "evidence_id": evidence_id, "study_id": study_id, "feature_id_standardized": feature,
        "phenotype": phenotype, "feature_type": feature_type, "simulated": False,
        "species_taxid": 29159, "mapping_confidence": "exact", "quality_flags": [],
        "adjusted_p_value": adjusted_p,
    }


def test_multi_omics_requires_a_second_layer_significant_under_the_same_phenotype():
    evidence = pd.DataFrame([
        _record("rna", "RNA", "gene", "thermal_tolerance", 1e-4),
        _record("meth", "METH", "methylation_region", "thermal_tolerance", 0.01),
    ])
    candidate = build_candidates(_gene_candidate(), evidence).iloc[0]
    assert candidate["tier"] == "multi_omics_convergence"
    assert candidate["n_supporting_layers"] == 2
    assert candidate["supporting_layers"] == "dna_methylation|transcriptomics"


def test_other_layer_evidence_under_a_different_phenotype_is_not_convergence():
    """A gene expressed under heat and methylated under pathogen challenge is
    two observations about two questions."""
    evidence = pd.DataFrame([
        _record("rna", "RNA", "gene", "thermal_tolerance", 1e-4),
        _record("meth", "METH", "methylation_region", "disease_susceptibility", 0.001),
    ])
    candidate = build_candidates(_gene_candidate(), evidence).iloc[0]
    assert candidate["tier"] == "emerging"
    assert candidate["n_supporting_layers"] == 1


def test_other_layer_evidence_must_itself_be_significant():
    evidence = pd.DataFrame([
        _record("rna", "RNA", "gene", "thermal_tolerance", 1e-4),
        _record("meth", "METH", "methylation_region", "thermal_tolerance", 0.4),
    ])
    candidate = build_candidates(_gene_candidate(), evidence).iloc[0]
    assert candidate["tier"] == "emerging"


def test_two_views_of_one_layer_do_not_count_as_two_layers():
    """A DMR and a single-CpG DML for the same gene are one methylation layer."""
    meta = _gene_candidate(feature_type="methylation_region", contributing_evidence_ids="dmr")
    evidence = pd.DataFrame([
        _record("dmr", "METH", "methylation_region", "thermal_tolerance", 1e-4),
        _record("dml", "METH", "methylation_locus", "thermal_tolerance", 1e-3),
    ])
    candidate = build_candidates(meta, evidence).iloc[0]
    assert candidate["n_supporting_layers"] == 1
    assert candidate["tier"] == "emerging"


def test_convergence_requires_the_candidate_itself_to_be_significant():
    """Two other layers converging on a gene this candidate shows no signal for
    is not evidence for this candidate."""
    meta = _gene_candidate(p_value=0.6, adjusted_p_value=0.9)
    evidence = pd.DataFrame([
        _record("rna", "RNA", "gene", "thermal_tolerance", 0.9),
        _record("meth", "METH", "methylation_region", "thermal_tolerance", 0.01),
        _record("prot", "PROT", "protein", "thermal_tolerance", 0.01),
    ])
    candidate = build_candidates(meta, evidence).iloc[0]
    assert candidate["n_supporting_layers"] == 2
    assert candidate["tier"] == "emerging"


def test_every_feature_type_declares_a_molecular_layer():
    from common import load_vocab

    for term_id, term in load_vocab("feature_types").items():
        assert term.get("molecular_layer"), f"{term_id} has no molecular_layer"


def test_demo_cross_layer_overlaps_are_all_cross_phenotype(isolated_reports):
    """The simulated studies share genes across RNA-seq, methylation and
    proteomics, but never under one phenotype, so none qualifies as
    convergence. The card must still show the adjacent evidence."""
    harmonize_all_demo_studies()
    meta_df = run_meta_analysis(feature_type="gene")
    candidates = build_candidates(meta_df, load_evidence_table())
    assert (candidates["tier"] != "multi_omics_convergence").all()
    assert (candidates["n_supporting_layers"] <= 1).all()


def test_score_does_not_saturate_above_the_old_ceilings():
    """Two candidates that both clear the old caps must still be distinguishable.

    Before 2026-09-05, `effect_magnitude_score` clipped at |log2FC| = 2 and
    `significance_score` at adjusted p = 1e-5, so 126 of the 1,994 high-priority
    candidates in the first real OsHV-1 pool scored an identical 68.57 — ties in
    exactly the region a reader picks validation targets from.
    """
    at_old_cap = _two_study_meta_row(pooled_effect=2.0, adjusted_p_value=1e-5)
    far_beyond = _two_study_meta_row(pooled_effect=5.8, adjusted_p_value=1e-75)

    score_at_cap = candidate_score(compute_components(at_old_cap, ["exact"], []))
    score_beyond = candidate_score(compute_components(far_beyond, ["exact"], []))

    assert score_beyond > score_at_cap
    # And the two components that used to clip are each strictly ordered.
    assert (compute_components(far_beyond, ["exact"], [])["effect_magnitude_score"]
            > compute_components(at_old_cap, ["exact"], [])["effect_magnitude_score"])
    assert (compute_components(far_beyond, ["exact"], [])["significance_score"]
            > compute_components(at_old_cap, ["exact"], [])["significance_score"])


def test_saturating_components_stay_within_bounds_at_extremes():
    extreme = _two_study_meta_row(
        k_studies=500, total_sample_size=10**6, pooled_effect=1e6, adjusted_p_value=1e-300,
        distinct_tissues=50, distinct_life_stages=50,
    )
    for key, value in compute_components(extreme, ["exact"], []).items():
        assert 0.0 <= value <= 1.0, f"{key}={value} out of bounds"


def test_disease_phenotype_is_not_labeled_resilience_evidence():
    """A strong tier must not be readable as a resilience claim.

    Both real OsHV-1 studies measure `disease_susceptibility` — infection
    response, with no survival or viral load recorded per animal — so their
    candidates are disease-associated evidence however well replicated.
    """
    candidate = build_candidates(
        pd.DataFrame([_two_study_meta_row(adjusted_p_value=1e-20, pooled_effect=3.0)]),
        _two_study_evidence(),
    ).iloc[0]

    assert candidate["tier"] == "high_priority_cross_study"
    assert candidate["evidence_class"] == "disease_associated"


def test_resilience_phenotype_is_labeled_resilience_evidence():
    candidate = build_candidates(
        pd.DataFrame([_two_study_meta_row(phenotype="survival", adjusted_p_value=1e-20)]),
        _two_study_evidence().assign(phenotype="survival"),
    ).iloc[0]
    assert candidate["evidence_class"] == "resilience_associated"


def test_context_replication_separates_study_count_from_context_breadth():
    evidence = _two_study_evidence()

    narrow = build_candidates(pd.DataFrame([_two_study_meta_row()]), evidence).iloc[0]
    assert narrow["k_studies"] == 2
    assert narrow["context_replication"] == "single_context"
    assert narrow["pooled_stressors"] == "pathogen_challenge"

    broad = build_candidates(
        pd.DataFrame([_two_study_meta_row(distinct_life_stages=2)]), evidence
    ).iloc[0]
    assert broad["context_replication"] == "multi_context"

    two_stressors = build_candidates(
        pd.DataFrame([_two_study_meta_row(distinct_stressors="pathogen_challenge|temperature")]),
        evidence,
    ).iloc[0]
    assert two_stressors["context_replication"] == "multi_context"

    single = build_candidates(
        pd.DataFrame([_two_study_meta_row(k_studies=1)]), evidence
    ).iloc[0]
    assert single["context_replication"] == "single_study"
