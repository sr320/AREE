"""Tests for orthogroup runs (src/mappings/orthogroups.py), WP6.

OrthoFinder itself is not run here; its output files are reproduced as
fixtures with the structural features that matter: a multi-copy group, a
single-copy group missing from one species, and unassigned genes.
"""
from __future__ import annotations

import json

import pytest

import mappings.orthogroups as og

SLUGS = {"mgigas": 29159, "mgall": 29158, "rphil": 129788}
ORTHOGROUPS = (
    "Orthogroup\tmgall\tmgigas\trphil\n"
    "OG0000000\t2001, 2002\t1001\t3001\n"      # expanded in mussel
    "OG0000001\t2003\t1002\t3002\n"            # one copy everywhere
    "OG0000002\t2004\t1003\t\n"                # absent from clam
)
UNASSIGNED = "Orthogroup\tmgall\tmgigas\trphil\nOG0000003\t\t1004\t\n"


@pytest.fixture
def results(tmp_path):
    d = tmp_path / "Results_Oct13" / "Orthogroups"
    d.mkdir(parents=True)
    (d / "Orthogroups.tsv").write_text(ORTHOGROUPS)
    (d / "Orthogroups_UnassignedGenes.tsv").write_text(UNASSIGNED)
    return d.parent


def test_parse_includes_unassigned_genes(results):
    rows = og.parse_orthogroups(results, SLUGS)
    assert ("2001", 29158, "OG0000000") in rows and ("2002", 29158, "OG0000000") in rows
    assert ("1004", 29159, "OG0000003") in rows
    assert len(rows) == 10


def test_a_gene_in_two_orthogroups_is_an_error(results):
    (results / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv").write_text(
        "Orthogroup\tmgall\tmgigas\trphil\nOG0000009\t\t1001\t\n")
    with pytest.raises(og.OrthogroupError, match="1001"):
        og.parse_orthogroups(results, SLUGS)


def test_an_unknown_species_column_is_an_error(results):
    with pytest.raises(og.OrthogroupError, match="not in this run"):
        og.parse_orthogroups(results, {"mgigas": 29159})


@pytest.mark.parametrize("copies,expected", [
    ({1: 1, 2: 1, 3: 1}, "single_copy_all"),
    ({1: 1, 2: 1, 3: 0}, "single_copy"),
    ({1: 2, 2: 1, 3: 1}, "multi_copy"),
    ({1: 3, 2: 0, 3: 0}, "species_specific"),
])
def test_copy_class(copies, expected):
    assert og.copy_class(copies, 3) == expected


def test_summary_and_written_run(results, tmp_path):
    rows = og.parse_orthogroups(results, SLUGS)
    taxids = [29158, 29159, 129788]
    summary = {r["orthogroup_id"]: r for r in og.summarize(rows, taxids)}
    assert summary["OG0000000"]["copy_class"] == "multi_copy"
    assert summary["OG0000000"]["copies_29158"] == 2
    assert summary["OG0000001"]["copy_class"] == "single_copy_all"
    assert summary["OG0000002"]["copy_class"] == "single_copy"
    assert summary["OG0000003"]["copy_class"] == "species_specific"

    staged = [{"slug": s, "ncbi_taxid": t, "proteome_sha256": "x" * 64, "proteins": 4}
              for s, t in sorted(SLUGS.items(), key=lambda kv: kv[1])]
    prov = og.write_run(tmp_path / "run1", rows, og.summarize(rows, [s["ncbi_taxid"] for s in staged]),
                        staged, version="OrthoFinder 3.1.5", command=["orthofinder", "-f", "x"],
                        results_dir=results, today="2026-10-13")
    on_disk = json.loads((tmp_path / "run1" / "provenance.json").read_text())
    assert on_disk == prov
    assert prov["orthogroups"] == {"total": 4, "by_copy_class": {
        "multi_copy": 1, "single_copy": 1, "single_copy_all": 1, "species_specific": 1}}
    assert prov["per_species"]["mgigas"] == {"genes": 4, "in_multi_species_orthogroups": 3}
    lines = (tmp_path / "run1" / "gene_orthogroups.tsv").read_text().splitlines()
    assert lines[0] == "ncbi_gene_id\tncbi_taxid\torthogroup_id" and len(lines) == 11


@pytest.fixture
def active_run(results, tmp_path, monkeypatch):
    """A written run in a temporary orthogroups directory, made active."""
    monkeypatch.setattr(og, "ORTHOGROUPS_DIR", tmp_path / "orthogroups")
    monkeypatch.setattr(og, "ACTIVE_PATH", tmp_path / "orthogroups" / "active.yaml")
    rows = og.parse_orthogroups(results, SLUGS)
    staged = [{"slug": s, "ncbi_taxid": t} for s, t in SLUGS.items()]
    og.write_run(tmp_path / "orthogroups" / "testrun", rows, og.summarize(rows, list(SLUGS.values())),
                 staged, version="OrthoFinder 3.1.5", command=None, results_dir=results)
    og.activate("testrun", tmp_path / "orthogroups" / "active.yaml")
    yield
    og.gene_orthogroup.cache_clear()


def test_lookup_through_the_active_run(active_run):
    assert og.gene_orthogroup("2002") == "OG0000000"
    assert og.gene_orthogroup(1002) == "OG0000001"
    assert og.gene_orthogroup("999999") is None
    assert og.active_run_ref()["run"] == "testrun"


def test_no_active_run_means_no_orthogroups(tmp_path, monkeypatch):
    monkeypatch.setattr(og, "ACTIVE_PATH", tmp_path / "absent.yaml")
    assert og.gene_orthogroup("1001") is None
    assert og.active_run_ref() is None


def test_resolution_attaches_the_orthogroup_for_a_real_crosswalk(active_run, tmp_path):
    from harmonize.identifiers import resolve_identifier, set_crosswalk_path

    crosswalk = tmp_path / "test_gene_id_crosswalk.tsv"
    crosswalk.write_text("ncbi_gene_id\tlocus_id\tgene_symbol\torthogroup_id\n2003\tLOC2003\tLOC2003\t\n")
    set_crosswalk_path(crosswalk)
    try:
        resolved = resolve_identifier("LOC2003", "locus_id")
        assert (resolved.feature_id_standardized, resolved.orthogroup_id) == ("2003", "OG0000001")
    finally:
        set_crosswalk_path(None)


def test_demo_crosswalk_keeps_its_own_orthogroups(active_run):
    from harmonize.identifiers import resolve_identifier, set_crosswalk_path

    set_crosswalk_path(None)
    resolved = resolve_identifier("LOC105333935", "ncbi_gene_id")
    assert resolved.orthogroup_id == "OG0000001"  # the demo file's illustrative value


def test_run_species_are_the_registered_reference_assemblies():
    species = {s["slug"]: s for s in og.run_species()}
    assert set(species) == {"mgigas", "mgall", "medul", "mtros", "mcali", "mmerc", "rphil"}
    assert species["mmerc"]["assembly_id"] == "MADL_Memer_1"  # the reference, not ASM1480567v1.1
    assert all(s["proteome"].name == "primary_proteins.faa" for s in species.values())


def test_staging_requires_each_proteome(tmp_path):
    missing = [{"slug": "x", "ncbi_taxid": 1, "ncbi_assembly_accession": "GCF_1.1",
                "proteome": tmp_path / "nope.faa"}]
    with pytest.raises(og.OrthogroupError, match="fetch-reference"):
        og.stage_proteomes(missing, tmp_path / "staged")


# --------------------------------------------------------------------------- #
# Invariants of the committed, active run
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def committed_run():
    run_dir = og.active_run_dir()
    if run_dir is None:
        pytest.skip("no active orthogroup run")
    return run_dir, json.loads((run_dir / "provenance.json").read_text())


def test_active_run_outputs_match_their_provenance(committed_run):
    from common import sha256_file

    run_dir, prov = committed_run
    for name, sha in prov["outputs"].items():
        assert sha256_file(run_dir / name) == sha, f"{name} changed since the run was written"


def test_active_run_covers_every_registered_species(committed_run):
    _, prov = committed_run
    assert {s["slug"] for s in prov["species"]} == {s["slug"] for s in og.run_species()}
    for sp in prov["species"]:
        assert prov["per_species"][sp["slug"]]["genes"] == sp["proteins"], sp["slug"]


def test_every_gene_in_the_active_run_is_in_its_species_crosswalk(committed_run):
    import pandas as pd

    from common import registered_crosswalks

    run_dir, prov = committed_run
    genes = pd.read_csv(run_dir / "gene_orthogroups.tsv", sep="\t", dtype=str)
    in_run = {int(s["ncbi_taxid"]) for s in prov["species"]}
    for entry in registered_crosswalks():
        if int(entry["ncbi_taxid"]) not in in_run:
            continue  # e.g. a GFF-route species not yet in an orthogroup run
        ids = set(pd.read_csv(entry["path"], sep="\t", comment="#", dtype=str,
                              usecols=["ncbi_gene_id"])["ncbi_gene_id"])
        mine = genes[genes["ncbi_taxid"] == str(entry["ncbi_taxid"])]
        assert len(mine) > 0, entry["path"].name
        assert set(mine["ncbi_gene_id"]) <= ids, entry["path"].name


def test_a_gene_belongs_to_one_orthogroup(committed_run):
    import pandas as pd

    run_dir, _ = committed_run
    genes = pd.read_csv(run_dir / "gene_orthogroups.tsv", sep="\t", dtype=str)
    assert genes["ncbi_gene_id"].is_unique
