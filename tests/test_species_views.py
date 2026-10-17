"""Per-species views for the dashboard and app, WP8 (src/reporting/species.py)."""
from __future__ import annotations

import pandas as pd
from conftest import harmonize_all_demo_studies

from common import load_vocab
from reporting.species import species_identity_basis, species_summary

OYSTER, MUSSEL = 29159, 29158


def _row(rows, taxid):
    return next(r for r in rows if r["ncbi_taxid"] == taxid)


def test_a_synonym_counts_with_its_species():
    """Crassostrea gigas and Magallana gigas are one species (taxid 29159)."""
    studies = [{"study_id": "A", "species": "Crassostrea gigas", "simulated": False},
               {"study_id": "B", "species": "Magallana gigas", "simulated": False},
               {"study_id": "C", "species": "Magallana gigas", "simulated": True}]
    oyster = _row(species_summary(studies), OYSTER)
    assert (oyster["n_studies_real"], oyster["n_studies_simulated"]) == (2, 1)
    assert oyster["studies_real"] == ["A", "B"]


def test_every_registered_species_appears_even_without_studies():
    rows = species_summary([])
    assert {r["ncbi_taxid"] for r in rows} == {t["ncbi_taxid"] for t in load_vocab("species").values()}
    assert all(r["n_studies_real"] == 0 for r in rows)


def test_crosswalk_route_and_size_are_reported():
    rows = {r["species"]: r for r in species_summary([])}
    assert rows["Magallana gigas"]["crosswalk_route"] == "ncbi_gene"
    assert rows["Mytilus coruscus"]["crosswalk_route"] == "gff_homology"
    assert rows["Mytilus coruscus"]["crosswalk_genes"] == 58540
    assert rows["Mytilus edulis"]["complex"] == "mytilus_edulis_complex"


def test_evidence_counts_are_split_by_origin():
    studies = [{"study_id": "A", "species": "Magallana gigas", "simulated": False}]
    evidence = pd.DataFrame({"study_id": ["A", "A", "D"], "simulated": ["False", "False", "True"],
                             "species_taxid": [OYSTER, OYSTER, MUSSEL]})
    rows = species_summary(studies, evidence)
    assert _row(rows, OYSTER)["n_evidence_records_real"] == 2
    assert _row(rows, OYSTER)["n_studies_real_harmonized"] == 1
    assert _row(rows, MUSSEL)["n_evidence_records_simulated"] == 1


def test_species_identity_basis():
    assert species_identity_basis({"species_identity": {"basis": "genotyped"}}) == "genotyped"
    assert species_identity_basis({}) is None


def test_dashboard_carries_species_rows_and_identity(isolated_reports, tmp_path, monkeypatch):
    import reporting.dashboard as dashboard

    harmonize_all_demo_studies()
    monkeypatch.setattr(dashboard, "EVIDENCE_TABLE", isolated_reports["evidence_table_path"])
    monkeypatch.setattr(dashboard, "META_ANALYSIS_DIR", isolated_reports["meta_dir"])
    monkeypatch.setattr(dashboard, "CANDIDATES_TABLE", tmp_path / "absent.tsv")
    data = dashboard.build_dashboard_data(tmp_path / "data.json")
    mussel = next(s for s in data["species"] if s["ncbi_taxid"] == MUSSEL)
    assert mussel["n_studies_simulated"] >= 1 and mussel["n_evidence_records_simulated"] > 0
    study = next(s for s in data["studies"] if s["study_id"] == "MYTGAL_HEAT01")
    assert study["species_identity_basis"] == "collection_site"
    assert study["species_complex"] == "mytilus_edulis_complex"
    assert study["species_canonical"] == "Mytilus galloprovincialis"
    hesser = next(s for s in data["studies"] if s["study_id"] == "HESSER2024_VCOR")
    assert hesser["species"] == "Crassostrea gigas" and hesser["species_canonical"] == "Magallana gigas"


def test_registry_index_records_species_identity_basis(isolated_registry):
    import csv

    from intake.registry import rebuild_registry

    rebuild_registry()
    with open(isolated_registry, newline="") as fh:
        rows = {r["study_id"]: r for r in csv.DictReader(fh)}
    assert rows["MYTGAL_HEAT01"]["species_identity_basis"] == "collection_site"
    assert rows["HESSER2024_VCOR"]["species_identity_basis"] == ""
