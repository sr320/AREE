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


def test_every_candidate_is_listed_but_only_signals_get_a_card(isolated_reports):
    """A genome-wide pool yields a candidate per gene. All of them belong in
    candidates.tsv; only those with a significant signal (pooled or in a
    contributing study) get a markdown card by default."""
    harmonize_all_demo_studies()
    evidence = load_evidence_table()
    source = evidence[
        (evidence["feature_id_standardized"] == "LOC105333935")
        & (evidence["phenotype"] == "larval_viability")
    ].iloc[0]

    rows = []
    for study_id, effect in [("NULL_A", 0.10), ("NULL_B", 0.12)]:
        row = source.copy()
        row["evidence_id"] = f"{study_id}-row"
        row["study_id"] = study_id
        row["comparison_id"] = f"{study_id}-comparison"
        row["feature_id_standardized"] = "GENE_NULL"
        row["feature_id_original"] = "GENE_NULL"
        row["effect_size"] = effect
        row["standard_error"] = 0.3
        row["p_value"] = 0.7
        row["adjusted_p_value"] = 0.95
        rows.append(row)
    for study_id, effect in [("SIG_A", 1.5), ("SIG_B", 1.7)]:
        row = source.copy()
        row["evidence_id"] = f"{study_id}-row"
        row["study_id"] = study_id
        row["comparison_id"] = f"{study_id}-comparison"
        row["feature_id_standardized"] = "GENE_SIG"
        row["feature_id_original"] = "GENE_SIG"
        row["effect_size"] = effect
        row["standard_error"] = 0.2
        row["p_value"] = 1e-9
        row["adjusted_p_value"] = 1e-7
        rows.append(row)
    path = isolated_reports["evidence_table_path"]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)

    index = build_evidence_cards(phenotype="larval_viability", feature_type="gene")
    cards_dir = isolated_reports["cards_dir"]
    candidates = pd.read_csv(cards_dir / "candidates.tsv", sep="\t")

    assert set(candidates["feature_id_standardized"]) == {"GENE_NULL", "GENE_SIG"}
    assert [e["feature_id_standardized"] for e in index] == ["GENE_SIG"]
    written = candidates.set_index("feature_id_standardized")["card_file"]
    assert Path(written["GENE_SIG"]).exists()
    assert pd.isna(written["GENE_NULL"])

    everything = build_evidence_cards(
        phenotype="larval_viability", feature_type="gene", all_cards=True
    )
    assert {e["feature_id_standardized"] for e in everything} == {"GENE_NULL", "GENE_SIG"}


def test_conflicting_candidate_keeps_its_card_via_study_level_significance(isolated_reports):
    """sod1 pools to nothing (p = 0.81) because two studies disagree, but one of
    them reported adjusted p = 0.03, so the conflict is surfaced on a card."""
    harmonize_all_demo_studies()
    index = build_evidence_cards(phenotype="larval_viability", feature_type="gene")
    sod1 = next(e for e in index if e["feature_id_standardized"] == "LOC105331241")
    assert sod1["adjusted_p_value"] > 0.05
    assert Path(sod1["card_file"]).exists()


def test_card_reports_the_adjusted_pooled_p_and_family_size(isolated_reports):
    harmonize_all_demo_studies()
    index = build_evidence_cards(phenotype="larval_viability", feature_type="gene")
    hsp70 = next(e for e in index if e["feature_id_standardized"] == "LOC105333935")
    content = Path(hsp70["card_file"]).read_text()
    assert "BH-adjusted across the" in content
    assert "| adj. p |" in content


def test_card_shows_other_layer_evidence_and_says_it_does_not_count(isolated_reports):
    harmonize_all_demo_studies()
    index = build_evidence_cards(phenotype="thermal_tolerance", feature_type="gene")
    hsp70 = next(e for e in index if e["feature_id_standardized"] == "LOC105333935")
    content = Path(hsp70["card_file"]).read_text()
    assert "own layer is `transcriptomics`" in content
    # methylation evidence for hsp70 exists, but under disease_susceptibility
    assert "dna_methylation" in content
    assert "| no | yes | no |" in content or "| no |" in content
    assert "**How the layers were linked:**" not in content


def test_intergenic_methylation_regions_stay_in_the_methylation_layer(isolated_reports):
    """An intergenic DMR is still a methylation region. `genomic_region` is the
    vocabulary's non-methylation feature type (QTL/locus); the gene-mapping
    outcome belongs in mapping_confidence, not the feature type."""
    harmonize_all_demo_studies()
    evidence = load_evidence_table()
    methylation = evidence[evidence["study_id"] == "GIGAS_PATH03"]
    assert set(methylation["feature_type"]) == {"methylation_region"}
    intergenic = methylation[methylation["feature_id_original"].str.startswith("DMR")]
    assert len(intergenic) == 2
    assert set(intergenic["mapping_confidence"]) == {"unresolved"}
    assert intergenic["feature_id_standardized"].isna().all()
    assert "genomic_region" not in set(evidence["feature_type"])
