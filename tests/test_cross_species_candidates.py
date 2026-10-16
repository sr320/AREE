"""Cross-species candidates and evidence cards, WP6 (src/prioritize/cross_species.py,
src/reporting/cross_species_cards.py)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from conftest import harmonize_all_demo_studies

from meta_analysis.orthogroups import run_orthogroup_meta_analysis
from prioritize.cross_species import (
    WEIGHTS,
    build_cross_species_candidates,
    compute_components,
    cross_species_score,
    species_estimates_for,
)
from reporting.cross_species_cards import (
    card_markdown,
    has_significant_signal,
    species_rows_for_card,
)

OYSTER, MUSSEL = 29159, 29158
CLASSES = {"OG1": "single_copy_all"}


def _gene(feature, taxid, effect, se, k, p, classification="resilience"):
    ids = "|".join(f"E{feature}_{i}" for i in range(k))
    return {"feature_id_standardized": feature, "phenotype": "thermal_tolerance", "feature_type": "gene",
            "simulated": False, "species_taxid": taxid, "k_studies": k,
            "studies": "|".join(f"S{feature}_{i}" for i in range(k)),
            "pooled_effect": effect, "pooled_se": se, "p_value": p, "adjusted_p_value": p,
            "contributing_evidence_ids": ids, "total_sample_size": 12 * k, "distinct_tissues": 1,
            "distinct_life_stages": 1, "distinct_stressors": "temperature",
            "mapping_confidences": "exact", "quality_flags_union": "", "_class": classification}


def _rank(*genes):
    genes = pd.DataFrame(genes)
    evidence = pd.DataFrame([
        {"evidence_id": eid, "feature_id_standardized": g["feature_id_standardized"],
         "feature_id_original": f"LOC{g['feature_id_standardized']}", "simulated": "False",
         "species_taxid": str(g["species_taxid"]), "orthogroup_id": "OG1",
         "resilience_classification": g["_class"], "study_id": eid.split("_")[0],
         "tissue": "gill", "life_stage": "adult", "stressor": "temperature"}
        for g in genes.to_dict("records") for eid in g["contributing_evidence_ids"].split("|")])
    result, species, _ = run_orthogroup_meta_analysis(genes.drop(columns="_class"), evidence, CLASSES)
    return build_cross_species_candidates(result, species, evidence), species, evidence


def test_two_replicated_concordant_species_converge():
    ranked, _, _ = _rank(_gene("101", OYSTER, 1.2, 0.2, k=2, p=1e-4), _gene("201", MUSSEL, 0.9, 0.25, k=2, p=1e-3))
    row = ranked.iloc[0]
    assert row["tier"] == "cross_species_convergence" and row["gates_failed"] == ""
    assert row["evidence_class"] == "resilience_associated"
    assert row["min_k_studies_per_species"] == 2 and row["total_sample_size"] == 48


def test_one_study_in_one_species_is_emerging_with_the_gate_named():
    ranked, _, _ = _rank(_gene("101", OYSTER, 1.2, 0.2, k=2, p=1e-4), _gene("201", MUSSEL, 0.9, 0.25, k=1, p=1e-3))
    row = ranked.iloc[0]
    assert row["tier"] == "cross_species_emerging" and row["gates_failed"] == "replication"


def test_exposure_only_evidence_in_one_species_blocks_convergence():
    """The weakest evidence class across all species' records decides, as at gene level."""
    ranked, _, _ = _rank(_gene("101", OYSTER, 1.2, 0.2, k=2, p=1e-4),
                         _gene("201", MUSSEL, 0.9, 0.25, k=2, p=1e-3, classification="exposure_only"))
    row = ranked.iloc[0]
    assert row["tier"] == "cross_species_emerging" and row["gates_failed"] == "phenotype_relevance"
    assert row["evidence_class"] == "exposure_only"


def test_opposite_signs_are_discordant_and_still_get_a_card():
    ranked, species, evidence = _rank(_gene("101", OYSTER, 1.5, 0.2, k=2, p=1e-6),
                                      _gene("201", MUSSEL, -1.5, 0.2, k=2, p=1e-6))
    row = ranked.iloc[0].to_dict()
    assert row["tier"] == "cross_species_discordant"
    assert "direction" in row["gates_failed"].split("|")
    rows = species_rows_for_card(row, species, evidence)
    # Pooled, the opposite effects cancel; the species' own significance keeps the card.
    assert row["adjusted_p_value"] > 0.05 and has_significant_signal(row, rows, 0.05)
    text = card_markdown(row, rows, "Thermal tolerance", {"run": "testrun"})
    assert "disagree in sign" in text and "Species disagree in direction" in text


def test_concordant_but_not_significant_is_unsupported():
    ranked, _, _ = _rank(_gene("101", OYSTER, 0.1, 0.5, k=2, p=0.8), _gene("201", MUSSEL, 0.05, 0.5, k=2, p=0.9))
    assert ranked.iloc[0]["tier"] == "cross_species_unsupported"


def test_score_is_reproducible_and_penalizes_between_species_heterogeneity():
    assert sum(WEIGHTS.values()) == 100
    agree, species, _ = _rank(_gene("101", OYSTER, 1.0, 0.1, k=2, p=1e-4), _gene("201", MUSSEL, 1.0, 0.1, k=2, p=1e-4))
    differ, _, _ = _rank(_gene("101", OYSTER, 0.4, 0.1, k=2, p=1e-4), _gene("201", MUSSEL, 1.6, 0.1, k=2, p=1e-4))
    row = agree.iloc[0].to_dict()
    components = compute_components(row, species_estimates_for(species, row), "resilience")
    assert cross_species_score(components) == row["score"]
    assert differ.iloc[0]["component_heterogeneity_penalty"] > 0.5
    assert differ.iloc[0]["score"] < agree.iloc[0]["score"]


def test_card_shows_each_species_side_by_side():
    ranked, species, evidence = _rank(_gene("101", OYSTER, 1.2, 0.2, k=2, p=1e-4),
                                      _gene("201", MUSSEL, 0.9, 0.25, k=1, p=1e-3))
    row = ranked.iloc[0].to_dict()
    text = card_markdown(row, species_rows_for_card(row, species, evidence), "Thermal tolerance",
                         {"run": "orthofinder_testrun"})
    assert "*Magallana gigas*" in text and "*Mytilus galloprovincialis*" in text
    assert "LOC101" in text and "LOC201" in text
    assert "orthofinder_testrun" in text
    assert "replication (fewer than two species have k >= 2" in text
    assert "not met" in text


def test_demo_builds_cross_species_cards_without_claiming_convergence(isolated_reports):
    from click.testing import CliRunner

    from aree.cli import main

    harmonize_all_demo_studies()
    out = CliRunner().invoke(main, ["build-cross-species-cards", "--feature-type", "gene"])
    assert out.exit_code == 0, out.output
    cards_dir = isolated_reports["cards_dir"] / "cross_species"
    ranked = pd.read_csv(cards_dir / "candidates.tsv", sep="\t")
    assert len(ranked) > 0
    # The demo mussel contributes one study per gene, so nothing converges.
    assert "cross_species_convergence" not in set(ranked["tier"])
    assert set(ranked["copy_class"]) == {"demo_crosswalk"}
    cards = sorted(cards_dir.glob("*.md"))
    assert cards and len(cards) == (ranked["card_file"].fillna("") != "").sum()
    for card in cards:
        text = card.read_text()
        assert "*Magallana gigas*" in text and "*Mytilus" in text
        assert "simulated demo evidence" in text and "not from an OrthoFinder run" in text


def test_demo_gene_cards_name_their_orthogroup(isolated_reports):
    from reporting.evidence_cards import build_evidence_cards_report

    harmonize_all_demo_studies()
    report = build_evidence_cards_report(phenotype="larval_viability", feature_type="gene")
    hsp70 = [r for r in report.index if r["feature_id_standardized"] == "LOC105333935"]
    assert hsp70, "demo HSP70 card missing"
    text = Path(hsp70[0]["card_file"]).read_text()
    assert "orthogroup `OG0000001`" in text


def test_dashboard_skips_orthogroup_results(isolated_reports, monkeypatch):
    """Regression: the gene-level pool summary globbed *_meta_analysis.tsv, which
    also matches the orthogroup output, whose columns differ."""
    import reporting.dashboard as dashboard
    from meta_analysis.run import write_meta_analysis, write_orthogroup_meta_analysis

    monkeypatch.setattr(dashboard, "META_ANALYSIS_DIR", isolated_reports["meta_dir"])
    harmonize_all_demo_studies()
    write_meta_analysis(None, "gene")
    write_orthogroup_meta_analysis(None, "gene")
    pools = dashboard._pool_summary()
    assert pools and all("k_studies" not in str(p.get("error", "")) for p in pools)


@pytest.mark.parametrize("copy_class,expected", [("single_copy_all", 1.0), ("multi_copy", 0.3)])
def test_orthology_score_reflects_copy_class(copy_class, expected):
    row = {"n_species": 2, "pooled_effect": 1.0, "p_value": 0.01, "adjusted_p_value": 0.01,
           "direction_consistency_species": 1.0, "copy_class": copy_class, "i_squared_between_species": 0.0}
    species = pd.DataFrame({"k_studies": [2, 2], "total_sample_size": [10, 10],
                            "mapping_confidences": ["exact", "exact"], "quality_flags_union": ["", ""]})
    assert compute_components(row, species, "resilience")["orthology_score"] == expected
