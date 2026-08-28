from pathlib import Path

from conftest import harmonize_all_demo_studies

from reporting.evidence_cards import build_evidence_cards


def test_build_evidence_cards_generates_files(isolated_reports):
    harmonize_all_demo_studies()
    index = build_evidence_cards(phenotype="larval_viability", feature_type="gene")
    assert len(index) > 0

    for entry in index:
        card_path = Path(entry["card_file"])
        assert card_path.exists()
        content = card_path.read_text()
        assert "Evidence card" in content
        assert "Recommended next validation step" in content
        assert "not validated biomarkers" in content


def test_evidence_card_flags_conflicting_direction(isolated_reports):
    harmonize_all_demo_studies()
    index = build_evidence_cards(phenotype="larval_viability", feature_type="gene")
    sod1_entry = next(e for e in index if e["feature_id_standardized"] == "LOC105331241")
    content = Path(sod1_entry["card_file"]).read_text()
    assert "Conflicting direction" in content


def test_every_candidate_gets_a_card_including_emerging(isolated_reports):
    harmonize_all_demo_studies()
    index = build_evidence_cards(phenotype="thermal_tolerance", feature_type="gene")
    tiers = {row["tier"] for row in index}
    assert "emerging" in tiers or "multi_omics_convergence" in tiers
