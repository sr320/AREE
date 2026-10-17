"""Curation tests for CRANDALL2022_HEAT_SEED.

Real proteomics study (PRIDE PXD015434; Crandall et al. 2022, PeerJ)
registered from PRIDE, the authors' Zenodo release and the hatchery notebook
in it. These tests pin the design the QC decision rests on (four pools: the
23 C control silo and the 29 C control silo, each on two days, four m/z-range
injections per pool), the sign of the authors' MSstats contrast, and the
decisions that follow: qc failed, nothing harmonized, biological_replicates 1.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "CRANDALL2022_HEAT_SEED.yaml"
STUDY_DIR = DATA_DIR / "studies" / "CRANDALL2022_HEAT_SEED"


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def samplesheet():
    with open(STUDY_DIR / "samplesheet.tsv") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def provenance():
    return json.loads((STUDY_DIR / "pride_provenance.json").read_text())


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_16_runs_are_4_pools_by_4_mz_range_injections(samplesheet, provenance):
    assert provenance["accession"] == "PXD015434"
    assert provenance["n_runs"] == len(samplesheet) == 16
    assert provenance["pride_files_by_category"] == {"RAW": 16, "SEARCH": 16}
    pools = Counter(r["pool_id"] for r in samplesheet)
    assert pools == dict.fromkeys(("Oysterseed_01", "Oysterseed_02", "Oysterseed_13", "Oysterseed_14"), 4)
    for pool in pools:
        groups = sorted(r["injection_group"] for r in samplesheet if r["pool_id"] == pool)
        assert groups == ["1", "2", "3", "4"]
    # The February 2016 re-acquisitions are one injection per pool, the same group.
    feb = {r["pool_id"]: r["injection_group"] for r in samplesheet if r["run_id"].startswith("2016_Feb")}
    assert set(feb) == set(pools) and len(set(feb.values())) == 1


def test_one_silo_per_arm_sampled_on_days_5_and_8(samplesheet, provenance):
    silos = {t: {r["silo"] for r in samplesheet if r["temperature_C"] == t} for t in ("23", "29")}
    assert silos == {"23": {"9"}, "29": {"8"}}
    days = Counter((r["temperature_C"], r["day_of_exposure"]) for r in samplesheet)
    assert set(days) == {(t, d) for t in ("23", "29") for d in ("5", "8")}
    assert {p["date"] for p in provenance["pools"].values()} == {"2015-09-11", "2015-09-14"}
    for r in samplesheet:
        assert r["condition"] == {"23": "control", "29": "treatment"}[r["temperature_C"]]
        assert r["msstats_bioreplicate"] == {"5": "1", "8": "2"}[r["day_of_exposure"]]


def test_msstats_contrast_is_29_minus_23_despite_its_label(provenance, study):
    ms = provenance["msstats_result"]
    assert ms["label"] == ["23C-29C"]
    assert ms["contrast_vector_23C_29C"] == [-1, 1]
    assert ms["n_proteins"] == 2808
    assert ms["n_adj_p_below_0.05"] == 0
    assert "reports 29C - 23C" in study["limitations"]


def test_failed_qc_nothing_harmonized(study):
    assert study["qc_status"] == "failed"
    assert study["analysis_status"] == "not_started"
    (comp,) = study["comparisons"]
    assert comp["results_file"] is None
    assert comp["biological_replicates"] == 1
    assert comp["sample_size"] == 4
    assert not list(STUDY_DIR.glob("*_standardized.tsv"))
    for flag in ("low_replication", "tank_level_replication", "batch_confounded"):
        assert flag in study["quality_flags"]


def test_survival_contrast_is_recorded_as_confounded(study):
    (comp,) = study["comparisons"]
    assert comp["comparison_id"] == "temp29_vs_temp23_days5and8"
    assert comp["meta_analysis_primary"] is True
    assert comp["phenotype"] == "survival"
    assert comp["control_condition"].startswith("23 C, Silo 9")
    assert comp["treatment_condition"].startswith("29 C, Silo 8")
    assert "ciliates" in comp["phenotype_unit"]
    assert "not evidence of a temperature effect on survival" in study["limitations"]


def test_zenodo_release_is_located_but_not_committed(provenance):
    assert provenance["zenodo"]["md5"] == "b5b97b6a8f0f805720263312e516e9d9"
    assert not list(STUDY_DIR.glob("*.zip")) and not list(STUDY_DIR.glob("*.xlsx"))


def test_samplesheet_checksum_matches_provenance(provenance):
    rel = "data/studies/CRANDALL2022_HEAT_SEED/samplesheet.tsv"
    digest = hashlib.sha256((REPO_ROOT / rel).read_bytes()).hexdigest()
    assert provenance["outputs"][rel] == digest
