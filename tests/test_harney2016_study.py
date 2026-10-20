"""Curation tests for HARNEY2016_OA_WARM_LARVAE.

Real proteomics study (PRIDE PXD002316; Harney et al. 2016, J Proteomics),
registered but deliberately not harmonized: a tank-replicated OA x warming
2-DE design whose spot volumes were never deposited, whose only published
results are 12 significance-selected spots without standard errors, and whose
UniProt accessions are mostly deleted. These tests pin those facts and the
decision, so a later curator does not harmonize Table 2 by accident.
"""
from __future__ import annotations

import json

import pytest
import yaml

from common import DATA_DIR, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "HARNEY2016_OA_WARM_LARVAE.yaml"
STUDY_DIR = DATA_DIR / "studies" / "HARNEY2016_OA_WARM_LARVAE"


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


def test_pride_holds_spot_identifications_only(provenance, study):
    assert provenance["accession"] == "PXD002316"
    assert provenance["quantitative_data_deposited"] is False
    assert provenance["n_pride_spot_runs"] == 22
    assert provenance["pride_files_by_category"] == {"OTHER": 22, "PEAK": 44, "RAW": 22, "RESULT": 22}
    assert study["data_availability"]["raw_data_available"] is False


def test_twelve_paper_spots_and_their_uniprot_status(provenance):
    spots = provenance["paper_spots"]
    assert sorted(spots, key=int) == [str(i) for i in range(1, 13)]
    assert spots["5"]["uniprot_accession"] == spots["6"]["uniprot_accession"] == "K1PE57"
    active = {s["uniprot_accession"]: s["ncbi_gene_ids"] for s in spots.values() if s["status"] == "active"}
    assert active == {"K1QWP8": ["105340468"], "K1Q350": ["105340512"]}
    assert provenance["uniprot_summary"] == {"deleted": 9, "active": 2}


def test_registered_not_harmonized(study):
    assert study["analysis_mode"] == "processed_results_harmonization"
    assert study["qc_status"] == "not_started"
    assert study["analysis_status"] == "not_started"
    assert all(c["results_file"] is None for c in study["comparisons"])
    assert not list(STUDY_DIR.glob("*.tsv"))
    assert not list(STUDY_DIR.glob("*.pdf"))
    for flag in ("processed_only", "no_multiple_testing_correction", "identifier_mapping_uncertain"):
        assert flag in study["quality_flags"]
    assert "not transcribed or harmonized" in study["limitations"]


def test_three_tanks_per_arm_and_five_table2_contrasts(study):
    comps = {c["comparison_id"]: c for c in study["comparisons"]}
    assert set(comps) == {
        "ph7_9_vs_ph8_1_at_20c", "ph7_9_vs_ph8_1_at_22c",
        "22c_vs_20c_at_ph8_1", "22c_vs_20c_at_ph7_9", "ph7_9_22c_vs_ph8_1_20c",
    }
    for c in comps.values():
        assert c["sample_size"] == 6 and c["biological_replicates"] == 3
    assert [cid for cid, c in comps.items() if c.get("meta_analysis_primary")] == ["ph7_9_vs_ph8_1_at_20c"]
    assert comps["ph7_9_22c_vs_ph8_1_20c"]["stressor_standardized"] == "multi_stressor"
    assert comps["22c_vs_20c_at_ph7_9"]["phenotype_direction"] == "not_applicable"
