"""Curation tests for CORPOREAU2014_OSHV_LOAD.

Real proteomics study (PRIDE PXD000905; Corporeau et al. 2014, J
Proteomics), registered but not harmonized: a 2-DE design (5 vs 5 oysters,
one gel each) whose spot volumes were never deposited; PRIDE holds only the
identification runs of the 25 excised spots. These tests pin those facts.
"""
from __future__ import annotations

import json

import pytest
import yaml

from common import DATA_DIR, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "CORPOREAU2014_OSHV_LOAD.yaml"
STUDY_DIR = DATA_DIR / "studies" / "CORPOREAU2014_OSHV_LOAD"


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def provenance():
    return json.loads((STUDY_DIR / "pride_provenance.json").read_text())


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_pride_holds_25_spot_identifications_only(provenance, study):
    assert provenance["quantitative_data_deposited"] is False
    assert provenance["n_spot_runs"] == 25
    assert provenance["pride_files_by_category"] == {"OTHER": 24, "PEAK": 24, "RAW": 25, "RESULT": 24, "SEARCH": 1}
    assert study["data_availability"]["raw_data_available"] is False


def test_spot_source_gels_are_the_table_1_oysters(provenance, study):
    sources = {tuple(x) for x in provenance["source_oysters_of_excised_spots"]}
    high = {("high_load_challenged", n) for n in (34, 37, 38, 42, 45)}
    low = {("low_load_control", n) for n in (61, 63, 67, 73, 79)}
    assert sources <= high | low
    (comp,) = study["comparisons"]
    for n in (34, 37, 38, 42, 45):
        assert str(n) in comp["treatment_condition"]
    for n in (61, 63, 67, 73, 79):
        assert str(n) in comp["control_condition"]


def test_registered_not_harmonized(study):
    assert study["qc_status"] == "not_started" and study["analysis_status"] == "not_started"
    (comp,) = study["comparisons"]
    assert comp["results_file"] is None
    assert (comp["sample_size"], comp["biological_replicates"]) == (10, 5)
    assert comp["resilience_classification"] == "disease"
    assert "tank_level_replication" in study["quality_flags"]
    assert not list(STUDY_DIR.glob("*.pdf"))
