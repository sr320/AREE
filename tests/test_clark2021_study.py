"""Curation tests for CLARK2021_HEAT_GREENLAND.

The first real mussel study: Greenland Mytilus edulis gill, a 2 h heat shock at
32 C against in-situ subtidal controls (primary) and against 22 C. These tests
pin the design read from ENA and ArrayExpress, the species-identity rule for
the M. edulis complex, the prespecified comparisons, and the two runs that
produced the committed results.
"""
from __future__ import annotations

import csv
import json
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_ID = "CLARK2021_HEAT_GREENLAND"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
RUN = ["heat_32c_vs_subtidal_control", "heat_32c_vs_22c"]


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def samplesheet():
    with open(STUDY_DIR / "samplesheet.tsv", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_five_groups_one_library_per_mussel(samplesheet):
    assert Counter(r["condition"] for r in samplesheet) == {
        "collected_from_outer_fjord_none": 5,
        "collected_from_inner_fjord_on_warm_day_none": 5,
        "collected_from_outer_fjord_on_warm_day_none": 5,
        "collected_from_outer_fjord_warmed_to_22_degrees": 5,
        "collected_from_outer_fjord_warmed_to_32_degrees": 4,
    }
    assert len({r["sample_accession"] for r in samplesheet}) == 24
    assert {r["library_layout"] for r in samplesheet} == {"PAIRED"}


def test_species_identity_is_declared_and_flagged(study):
    """M. edulis is in the hybridizing complex; collection site is not genotyping."""
    assert study["species_identity"]["basis"] == "collection_site"
    assert "species_assignment_unverified" in study["quality_flags"]


def test_one_primary_exposure_comparison(study):
    comps = {c["comparison_id"]: c for c in study["comparisons"]}
    assert [cid for cid, c in comps.items() if c.get("meta_analysis_primary")] == ["heat_32c_vs_subtidal_control"]
    for c in comps.values():
        assert c["phenotype"] == "thermal_tolerance" and c["tissue"] == "gill"
        assert c["resilience_classification"] == "stress_response"
        assert c["contrast_type"] == "exposure_response"


@pytest.mark.parametrize("cid", RUN)
def test_each_run_is_the_one_registered(study, cid):
    comp = {c["comparison_id"]: c for c in study["comparisons"]}[cid]
    assert (REPO_ROOT / comp["results_file"]).exists()
    manifest = json.loads((STUDY_DIR / f"{STUDY_ID}_{cid}_workflow_manifest.json").read_text())
    assert manifest["parameters"]["tag_seq"] is False
    assert manifest["deseq2_design"]["n_libraries_or_units"] == 9
    assert manifest["deseq2_design"]["residual_df"] == 7
    assert all(f.get("matches_bundle_provenance", True) for f in manifest["reference"]["files"])
    assert manifest["reference"]["genome_assembly"] == "GCF_963676685.1"
    assert manifest["warnings"] == []


def test_replicate_dispersion_qc_passes():
    with open(STUDY_DIR / "replicate_dispersion_qc.tsv", encoding="utf-8") as fh:
        rows = {r["condition"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert set(rows) == {"collected_from_outer_fjord_none", "collected_from_outer_fjord_warmed_to_22_degrees",
                         "collected_from_outer_fjord_warmed_to_32_degrees"}
    assert all(r["passes"] == "True" for r in rows.values())


def test_the_unrun_comparison_has_no_results_yet(study):
    comp = {c["comparison_id"]: c for c in study["comparisons"]}["inner_fjord_warm_day_vs_subtidal_control"]
    assert not (REPO_ROOT / comp["results_file"]).exists()
    assert study["analysis_status"] == "in_progress"
