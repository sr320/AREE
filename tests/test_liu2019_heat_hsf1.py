"""Curation tests for LIU2019_HEAT_HSF1.

Liu et al. 2019 (PRJNA516762): 35 C heat shock vs control in juvenile gill,
three pooled libraries (five oysters each) per group. Only the direct
heat-vs-control contrast is used; the HSF1-RNAi groups are not. These tests
pin the design, the single prespecified comparison, and the run outputs.
"""
from __future__ import annotations

import csv
import json

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "LIU2019_HEAT_HSF1.yaml"
STUDY_DIR = DATA_DIR / "studies" / "LIU2019_HEAT_HSF1"
CID = "heat_35c_2h_vs_control"


def _tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def samplesheet():
    return _tsv(STUDY_DIR / "samplesheet.tsv")


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_design_is_four_groups_of_three(samplesheet):
    prov = json.loads((STUDY_DIR / "ena_provenance.json").read_text())
    assert prov["bioproject"] == "PRJNA516762"
    assert prov["group_sizes"] == {"HS": 3, "HS_RNAi": 3, "RNAi": 3, "control": 3}
    for row in samplesheet:
        assert row["sample_alias"].startswith(row["treatment"] + "-")


def test_only_heat_vs_control_is_used(study):
    (comp,) = study["comparisons"]
    assert comp["comparison_id"] == CID
    assert comp["meta_analysis_primary"] is True
    assert comp["stressor_standardized"] == "temperature"
    assert comp["resilience_classification"] == "stress_response"
    assert comp["sample_size"] == 6 and comp["biological_replicates"] == 3


def test_run_outputs_match_the_registered_comparison(study):
    (comp,) = study["comparisons"]
    params = json.loads((STUDY_DIR / f"LIU2019_HEAT_HSF1_{CID}_workflow_manifest.json").read_text())["parameters"]
    assert params["comparison_id"] == CID
    assert params["control_levels"] == ["control"] and params["treatment_levels"] == ["HS"]
    rows = _tsv(REPO_ROOT / comp["results_file"])
    assert len(rows) > 20000
    assert sum(r["lfcSE"] not in ("", "NA") for r in rows) > 0.95 * len(rows)
    dispersion = _tsv(STUDY_DIR / "replicate_dispersion_qc.tsv")
    assert {r["condition"] for r in dispersion} == {"control", "HS"}
