from pathlib import Path

import pandas as pd
from conftest import harmonize_all_demo_studies

from harmonize.core import load_evidence_table
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


def test_card_files_and_evidence_are_partitioned_by_origin_and_species(isolated_reports):
    harmonize_all_demo_studies()
    evidence = load_evidence_table()
    source = evidence[
        (evidence["feature_id_standardized"] == "LOC105331241")
        & (evidence["phenotype"] == "larval_viability")
    ].iloc[0]

    rows = []
    for evidence_id, study_id, simulated, species_taxid in [
        ("simulated-row", "SIMULATED_STUDY", True, 29159),
        ("real-row", "REAL_STUDY", False, 29159),
        ("other-species-row", "OTHER_SPECIES", False, 99999),
    ]:
        row = source.copy()
        row["evidence_id"] = evidence_id
        row["study_id"] = study_id
        row["comparison_id"] = f"{evidence_id}-comparison"
        row["simulated"] = simulated
        row["species_taxid"] = species_taxid
        row["standard_error"] = 0.2
        row["p_value"] = 0.01
        rows.append(row)

    path = isolated_reports["evidence_table_path"]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)

    index = build_evidence_cards(phenotype="larval_viability", feature_type="gene")
    assert len(index) == 3
    assert len({entry["card_file"] for entry in index}) == 3
    for entry in index:
        content = Path(entry["card_file"]).read_text()
        expected_study = {
            (True, 29159): "SIMULATED_STUDY",
            (False, 29159): "REAL_STUDY",
            (False, 99999): "OTHER_SPECIES",
        }[(bool(entry["simulated"]), int(entry["species_taxid"]))]
        assert expected_study in content
        assert sum(name in content for name in [
            "SIMULATED_STUDY", "REAL_STUDY", "OTHER_SPECIES"
        ]) == 1
