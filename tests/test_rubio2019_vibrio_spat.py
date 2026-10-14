"""Curation tests for RUBIO2019_VIBRIO_SPAT.

Rubio et al. 2019 (PRJNA515169): spat injected with virulent or non-virulent
vibrios of two species, or mock, sampled at 8 h; one pool of 30 spat per
condition in each of three independent experiments. These tests pin the
design read from the aliases, the experiment blocking, and the four
prespecified comparisons.
"""
from __future__ import annotations

import csv
import json
import re

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "RUBIO2019_VIBRIO_SPAT.yaml"
STUDY_DIR = DATA_DIR / "studies" / "RUBIO2019_VIBRIO_SPAT"
CONFIG = REPO_ROOT / "config" / "RUBIO2019_VIBRIO_SPAT.config"


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


def test_six_conditions_by_three_experiments(samplesheet):
    prov = json.loads((STUDY_DIR / "ena_provenance.json").read_text())
    assert prov["bioproject"] == "PRJNA515169" and prov["n_runs"] == 18
    assert set(prov["group_sizes"].values()) == {3} and len(prov["group_sizes"]) == 6
    cells = {(r["condition"], r["experiment"]) for r in samplesheet}
    assert len(cells) == 18  # every condition in every experiment exactly once


def test_comparisons_are_blocked_on_experiment(study):
    assert re.search(r"covariates\s*=\s*'experiment'", CONFIG.read_text())
    ids = [c["comparison_id"] for c in study["comparisons"]]
    assert ids == ["vcrass_j29_virulent_vs_mock_8h", "vtasma_lgp32_virulent_vs_mock_8h",
                   "vcrass_virulent_vs_nonvirulent_8h", "vtasma_virulent_vs_nonvirulent_8h"]
    primary = [c["comparison_id"] for c in study["comparisons"] if c.get("meta_analysis_primary")]
    assert primary == ["vcrass_j29_virulent_vs_mock_8h"]


def test_pathogen_challenge_is_disease_evidence(study):
    for comp in study["comparisons"]:
        assert comp["stressor_standardized"] == "pathogen_challenge"
        assert comp["phenotype"] == "disease_susceptibility"
        assert comp["resilience_classification"] == "disease"
        assert comp["sample_size"] == 6 and comp["biological_replicates"] == 3


def test_run_outputs_match_the_registered_comparisons(study):
    expected = {
        "vcrass_j29_virulent_vs_mock_8h": ("anesthesis", "Vcrass_J2_9"),
        "vtasma_lgp32_virulent_vs_mock_8h": ("anesthesis", "Vtasma_LGP32"),
        "vcrass_virulent_vs_nonvirulent_8h": ("Vcrass_J2_8", "Vcrass_J2_9"),
        "vtasma_virulent_vs_nonvirulent_8h": ("Vtasma_LMG20012T", "Vtasma_LGP32"),
    }
    for comp in study["comparisons"]:
        cid = comp["comparison_id"]
        params = json.loads((STUDY_DIR / f"RUBIO2019_VIBRIO_SPAT_{cid}_workflow_manifest.json").read_text())["parameters"]
        assert params["comparison_id"] == cid
        assert (params["control_levels"], params["treatment_levels"]) == ([expected[cid][0]], [expected[cid][1]])
        assert params["covariates"] == ["experiment"]
        rows = _tsv(REPO_ROOT / comp["results_file"])
        assert len(rows) > 20000
        assert sum(r["lfcSE"] not in ("", "NA") for r in rows) > 0.95 * len(rows)
