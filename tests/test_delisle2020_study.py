"""Curation tests for DELISLE2020_OSHV_TEMP.

This study is AREE's first independent raw-reanalysis partner selected to make
real cross-study pathogen-challenge pooling possible once both raw workflows
have been run and harmonized.
"""
from __future__ import annotations

import csv
import json

import pandas as pd
import pytest
import yaml

from common import DATA_DIR, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "DELISLE2020_OSHV_TEMP.yaml"
STUDY_DIR = DATA_DIR / "studies" / "DELISLE2020_OSHV_TEMP"


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def samplesheet():
    with open(STUDY_DIR / "samplesheet.tsv") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def ena_provenance():
    return json.loads((STUDY_DIR / "ena_provenance.json").read_text())


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_samplesheet_matches_live_ena_design_snapshot(samplesheet, ena_provenance):
    assert ena_provenance["bioproject"] == "PRJNA593309"
    assert len(samplesheet) == 43
    assert ena_provenance["n_runs"] == 43
    assert ena_provenance["n_samples"] == 43


def test_registered_comparison_is_the_pool_compatible_pathogen_endpoint(study):
    assert [c["comparison_id"] for c in study["comparisons"]] == ["oshv1_21c_96h_vs_21c_0h"]
    comp = study["comparisons"][0]
    assert study["analysis_status"] == "harmonized"
    assert study["qc_status"] == "in_progress"
    assert comp["stressor_standardized"] == "pathogen_challenge"
    assert comp["phenotype"] == "disease_susceptibility"
    assert comp["resilience_classification"] == "disease"
    assert comp["sample_size"] == 6
    assert comp["biological_replicates"] == 3
    assert comp["results_file"] == (
        "data/studies/DELISLE2020_OSHV_TEMP/"
        "DELISLE2020_OSHV_TEMP_oshv1_21c_96h_vs_21c_0h_dge_standardized.tsv"
    )
    assert (DATA_DIR.parent / comp["results_file"]).exists()


def test_registered_comparison_groups_are_replicated(samplesheet):
    by_condition: dict[str, list[dict]] = {}
    for row in samplesheet:
        by_condition.setdefault(row["condition"], []).append(row)

    assert {len(by_condition["21_0"]), len(by_condition["21_96h"])} == {3}
    assert {row["temp"] for row in by_condition["21_0"] + by_condition["21_96h"]} == {"21"}
    assert {row["time_post_infection"] for row in by_condition["21_0"]} == {"0"}
    assert {row["time_post_infection"] for row in by_condition["21_96h"]} == {"96h"}


def test_low_replication_26c_timepoints_are_not_registered(study, ena_provenance):
    assert ena_provenance["group_sizes"]["26_12h"] == 2
    assert ena_provenance["group_sizes"]["26_24h"] == 2
    registered = " ".join(c["comparison_id"] for c in study["comparisons"])
    assert "26c" not in registered


def test_fastq_manifest_covers_every_run_as_paired_end(samplesheet):
    with open(STUDY_DIR / "fastq_manifest.tsv") as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    runs_in_sheet = {r["run_accession"] for r in samplesheet}
    runs_in_manifest = {r["run_accession"] for r in manifest}
    assert runs_in_sheet == runs_in_manifest
    assert len(manifest) == 2 * len(runs_in_sheet)
    assert all(len(r["md5"]) == 32 for r in manifest)


def test_committed_outputs_support_the_first_real_cross_study_pool(study):
    calla_path = DATA_DIR / "studies" / "CALLA2026_OSHV" / (
        "CALLA2026_OSHV_miyagi_oshv1_usa_vs_control_dge_standardized.tsv"
    )
    delisle_path = DATA_DIR / "studies" / "DELISLE2020_OSHV_TEMP" / (
        "DELISLE2020_OSHV_TEMP_oshv1_21c_96h_vs_21c_0h_dge_standardized.tsv"
    )
    calla = pd.read_csv(calla_path, sep="\t", usecols=["gene_id"], dtype=str)
    delisle = pd.read_csv(delisle_path, sep="\t", usecols=["gene_id"], dtype=str)

    shared_gene_ids = set(calla["gene_id"]) & set(delisle["gene_id"])
    assert len(shared_gene_ids) == 23094
    assert study["comparisons"][0]["phenotype"] == "disease_susceptibility"


def test_bioproject_mismatch_is_caught(tmp_path, monkeypatch):
    import intake.schema_validate as sv

    doc = yaml.safe_load(STUDY_PATH.read_text())
    doc["accessions"]["bioproject"] = "PRJNA9999999"
    path = tmp_path / "DELISLE2020_OSHV_TEMP.yaml"
    path.write_text(yaml.safe_dump(doc))
    monkeypatch.setattr(sv, "DATA_DIR", DATA_DIR)

    result = sv.validate_study_file(path)
    assert not result.valid
    assert any("PRJNA9999999" in e for e in result.errors)
