"""Curation tests for JIN2025_VIBRIO_MUCUS.

Real, unpublished proteomics dataset (PRIDE PXD064651), registered but not
reanalyzable: an isobaric (IBT 16-plex) experiment whose deposit names six
samples but no reporter-channel map, and whose Mascot search carried no
quantitation. These tests pin those facts and the register-only decision.
"""
from __future__ import annotations

import json

import pytest
import yaml

from common import DATA_DIR, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "JIN2025_VIBRIO_MUCUS.yaml"
STUDY_DIR = DATA_DIR / "studies" / "JIN2025_VIBRIO_MUCUS"


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


def test_isobaric_deposit_without_channel_map_or_quantification(provenance):
    assert provenance["accession"] == "PXD064651"
    assert provenance["submission_type"] == "PARTIAL"
    assert provenance["references"] == []
    assert provenance["channel_to_sample_key_deposited"] is False
    assert provenance["quantification_table_deposited"] is False
    assert provenance["mascot_dat_parameters"]["QUANTITATION"] == "None"
    assert "IBT16plex" in provenance["mascot_dat_parameters"]["MODS"]
    assert provenance["pride_files_by_category"] == {"PEAK": 1, "RAW": 3, "SEARCH": 1}


def test_registered_not_reanalyzed(study):
    assert study["qc_status"] == "not_started" and study["analysis_status"] == "not_started"
    (comp,) = study["comparisons"]
    assert comp["results_file"] is None
    assert comp["stressor_standardized"] == "pathogen_challenge"
    assert comp["resilience_classification"] == "exposure_only"
    assert (comp["sample_size"], comp["biological_replicates"]) == (6, 3)
    assert study["accessions"]["doi"] is None
    assert "No channel-to-sample key" in study["limitations"]
