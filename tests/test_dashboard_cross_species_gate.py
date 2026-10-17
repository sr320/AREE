"""The WP6 dashboard gate (src/reporting/dashboard.py).

No cross-species pooling is shown on the dashboard until at least two species
each have a within-species pool with k >= 2 on the same phenotype. Before that
the orthogroup view is a lookup, not a meta-analysis.
"""
from __future__ import annotations

import pandas as pd
from click.testing import CliRunner
from conftest import ALL_DEMO_STUDY_IDS

import reporting.dashboard as dashboard
from aree.cli import main
from reporting.dashboard import build_dashboard_data, cross_species_gate

OYSTER, MUSSEL, CLAM = 29159, 29158, 129788


def _pool(phenotype, taxid, max_k, simulated=False, feature_type="gene"):
    return {"phenotype": phenotype, "feature_type": feature_type, "simulated": simulated,
            "species_taxid": taxid, "species": dashboard._species_label(taxid), "max_k": max_k,
            "studies": [f"S{taxid}_{i}" for i in range(max_k)]}


def test_gate_opens_when_two_species_each_have_k2():
    (g,) = cross_species_gate([_pool("thermal_tolerance", OYSTER, 3), _pool("thermal_tolerance", CLAM, 2)])
    assert g["open"] is True and g["n_species_k2"] == 2
    assert [s["species"] for s in g["species"]] == ["Magallana gigas", "Ruditapes philippinarum"]


def test_gate_stays_closed_when_one_species_has_a_single_study():
    (g,) = cross_species_gate([_pool("thermal_tolerance", OYSTER, 3), _pool("thermal_tolerance", CLAM, 1)])
    assert g["open"] is False and "only 1 of 2 species" in g["reason"]


def test_gate_is_per_phenotype_and_per_origin():
    gate = cross_species_gate([
        _pool("thermal_tolerance", OYSTER, 2), _pool("thermal_tolerance", CLAM, 2),
        _pool("larval_viability", OYSTER, 2, simulated=True), _pool("larval_viability", MUSSEL, 1, simulated=True),
        _pool("disease_susceptibility", OYSTER, 4),          # one species: not a cross-species phenotype
        _pool("thermal_tolerance", MUSSEL, 5, feature_type="protein"),  # only gene pools count
    ])
    assert {(g["phenotype"], g["simulated"]): g["open"] for g in gate} == {
        ("thermal_tolerance", False): True, ("larval_viability", True): False}


def _write_cross_table(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_rows_for_a_closed_phenotype_are_withheld_not_ranked(tmp_path, monkeypatch):
    table = tmp_path / "cross_species" / "candidates.tsv"
    base = {"feature_type": "gene", "copy_class": "single_copy_all", "evidence_class": "exposure_only",
            "species_taxids": f"{OYSTER}|{CLAM}", "genes": f"{OYSTER}:1|{CLAM}:2",
            "k_studies_by_species": f"{OYSTER}:3|{CLAM}:2", "pooled_effect": 0.5, "adjusted_p_value": 0.01,
            "i_squared_between_species": 0.0, "simulated": False, "card_file": ""}
    _write_cross_table(table, [
        base | {"orthogroup_id": "OG1", "phenotype": "thermal_tolerance", "tier": "cross_species_emerging", "score": 50},
        base | {"orthogroup_id": "OG2", "phenotype": "acidification_tolerance", "tier": "cross_species_emerging",
                "score": 60},
    ])
    monkeypatch.setattr(dashboard, "CROSS_SPECIES_TABLE", table)
    gate = cross_species_gate([_pool("thermal_tolerance", OYSTER, 3), _pool("thermal_tolerance", CLAM, 2),
                               _pool("acidification_tolerance", OYSTER, 3),
                               _pool("acidification_tolerance", CLAM, 1)])
    summary = dashboard._cross_species_summary(gate, {})
    assert [r["orthogroup_id"] for r in summary["top_real"]] == ["OG1"]
    assert summary["withheld"] == [{"simulated": False, "phenotype": "acidification_tolerance", "n": 1}]
    (row,) = summary["top_real"]
    assert [(g["species"], g["gene"], g["k_studies"]) for g in row["species_genes"]] == [
        ("Magallana gigas", "1", 3), ("Ruditapes philippinarum", "2", 2)]


def test_demo_dashboard_gates_the_mussel_and_offers_the_lookup(isolated_reports, isolated_registry,
                                                                tmp_path, monkeypatch):
    """The demo mussel has one study, so its larval rows are withheld; the
    oyster candidates still show where their mussel orthologs were tested."""
    runner = CliRunner()
    for study_id in ALL_DEMO_STUDY_IDS:
        assert runner.invoke(main, ["register-study", f"registry/studies/{study_id}.yaml"]).exit_code == 0
        assert runner.invoke(main, ["harmonize", "--study", study_id, "--date", "2026-01-01"]).exit_code == 0
    assert runner.invoke(main, ["meta-analyze", "--feature-type", "gene"]).exit_code == 0
    assert runner.invoke(main, ["build-evidence-cards"]).exit_code == 0
    out = runner.invoke(main, ["build-cross-species-cards", "--feature-type", "gene"])
    assert out.exit_code == 0, out.output

    reports_dir = isolated_reports["evidence_table_path"].parent.parent
    monkeypatch.setattr(dashboard, "EVIDENCE_TABLE", reports_dir / "evidence" / "evidence_table.tsv")
    monkeypatch.setattr(dashboard, "CANDIDATES_TABLE", reports_dir / "evidence_cards" / "candidates.tsv")
    monkeypatch.setattr(dashboard, "CROSS_SPECIES_TABLE",
                        reports_dir / "evidence_cards" / "cross_species" / "candidates.tsv")
    monkeypatch.setattr(dashboard, "META_ANALYSIS_DIR", reports_dir / "meta_analysis")
    data = build_dashboard_data(tmp_path / "site" / "dashboard" / "data.json")

    larval = [g for g in data["cross_species"]["gate"] if g["phenotype"] == "larval_viability"]
    assert len(larval) == 1 and larval[0]["simulated"] is True and larval[0]["open"] is False
    assert data["cross_species"]["top_simulated"] == []
    assert data["cross_species"]["n_withheld"] > 0
    assert data["artifacts_present"]["cross_species_candidates"] is True

    # Pools are per species: the larval phenotype is two pools, not one.
    larval_pools = [p for p in data["pools"] if p["phenotype"] == "larval_viability"]
    assert {p["species_taxid"] for p in larval_pools} == {OYSTER, MUSSEL}

    rows = [r for rows in data["candidates"]["top_simulated_by_phenotype"].values() for r in rows]
    hsp70 = [r for r in rows if r["feature_id_standardized"] == "LOC105333935" and r["species_taxid"] == OYSTER]
    assert hsp70, "demo oyster HSP70 missing from the per-phenotype lists"
    (mussel,) = hsp70[0]["orthologs"]
    assert mussel["species_taxid"] == MUSSEL and "larval_viability" in mussel["phenotypes"]
    assert all(o["species_taxid"] != r["species_taxid"] for r in rows for o in r["orthologs"])
