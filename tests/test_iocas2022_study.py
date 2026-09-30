"""Curation tests for IOCAS2022_OA_ENERGY.

AREE's first real ocean-acidification study, registered from deposited
metadata before any reanalysis. These tests pin the design read from ENA and
the curation decisions that are easy to erode later: exposure_only
classification, one prespecified pooled comparison, and no borrowed DOI.
"""
from __future__ import annotations

import csv
import json

import pytest
import yaml

from common import DATA_DIR, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "IOCAS2022_OA_ENERGY.yaml"
STUDY_DIR = DATA_DIR / "studies" / "IOCAS2022_OA_ENERGY"
TIMEPOINTS = ("7d", "28d", "56d")


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
    assert ena_provenance["bioproject"] == "PRJNA826964"
    assert ena_provenance["n_runs"] == ena_provenance["n_samples"] == len(samplesheet) == 18
    assert ena_provenance["min_group_size"] == 3
    # One library per animal: this is what rejected PRJNA623063 lacked.
    assert len({r["sample_alias"] for r in samplesheet}) == 18


def test_each_comparison_is_against_its_time_matched_control(study, samplesheet):
    by_condition: dict[str, int] = {}
    for row in samplesheet:
        by_condition[row["condition"]] = by_condition.get(row["condition"], 0) + 1

    ids = [c["comparison_id"] for c in study["comparisons"]]
    assert ids == [f"oa_{t}_vs_control_{t}" for t in TIMEPOINTS]
    for t in TIMEPOINTS:
        assert by_condition[f"control_{t}"] == 3
        assert by_condition[f"ocean_acidification_{t}"] == 3


def test_no_measured_phenotype_means_exposure_only(study):
    for comp in study["comparisons"]:
        assert comp["stressor_standardized"] == "ocean_acidification"
        assert comp["resilience_classification"] == "exposure_only"
        assert comp["phenotype_direction"] == "not_applicable"
        assert comp["tissue"] == "digestive_gland"
        assert comp["sample_size"] == 6
        assert comp["results_file"] is None
    assert "ambiguous_phenotype_definition" in study["quality_flags"]


def test_exactly_one_prespecified_pooled_comparison(study):
    primary = [c["comparison_id"] for c in study["comparisons"] if c.get("meta_analysis_primary")]
    assert primary == ["oa_28d_vs_control_28d"]


def test_same_lab_gonad_paper_is_not_attached_as_the_source(study):
    assert study["accessions"]["doi"] is None
    assert "cirep" not in study["citation"]


def test_fastq_manifest_covers_every_run_as_paired_end(samplesheet):
    with open(STUDY_DIR / "fastq_manifest.tsv") as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    runs_in_sheet = {r["run_accession"] for r in samplesheet}
    assert {r["run_accession"] for r in manifest} == runs_in_sheet
    assert len(manifest) == 2 * len(runs_in_sheet)
    assert all(len(r["md5"]) == 32 for r in manifest)
