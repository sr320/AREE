"""Curation tests for WANAMAKER2020_HEAT_SEED.

Real proteomics study (PRIDE PXD013262; Wanamaker et al. 2020, BMC Genomics)
registered from PRIDE, the authors' figshare package and the paper. These
tests pin the design the QC decision rests on (one silo per temperature, one
pool per temperature per day, two technical MS runs per pool) and the
decisions that follow from it: qc failed, nothing harmonized, the phenotype
classified exposure_only, and biological_replicates = 1.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "WANAMAKER2020_HEAT_SEED.yaml"
STUDY_DIR = DATA_DIR / "studies" / "WANAMAKER2020_HEAT_SEED"
DAYS = (3, 5, 7, 9, 11, 13)


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


def test_26_runs_are_13_pools_in_technical_duplicate(samplesheet, provenance):
    assert provenance["accession"] == "PXD013262"
    assert provenance["n_runs"] == len(samplesheet) == 26
    pools = Counter(r["pool_id"] for r in samplesheet)
    assert len(pools) == 13 and set(pools.values()) == {2}
    for pool in pools:
        reps = sorted(r["technical_replicate"] for r in samplesheet if r["pool_id"] == pool)
        assert reps == ["1", "2"]
    assert provenance["pride_files_by_category"] == {"OTHER": 1, "PEAK": 26, "RAW": 26, "SEARCH": 26}


def test_one_silo_per_temperature_and_one_pool_per_day(samplesheet, provenance):
    silos = {t: {r["silo"] for r in samplesheet if r["temperature_C"] == t} for t in ("23", "29")}
    assert silos == {"23": {"3"}, "29": {"9"}}
    assert provenance["silos_by_temperature_C"] == {"23": ["3"], "29": ["9"]}
    cells = Counter((r["temperature_C"], r["day_of_exposure"]) for r in samplesheet if r["silo"])
    assert set(cells) == {(t, str(d)) for t in ("23", "29") for d in DAYS}
    assert set(cells.values()) == {2}  # two technical runs of one pool
    (baseline,) = {r["pool_id"] for r in samplesheet if r["condition"] == "baseline"}
    assert baseline == "pool0"


def test_dpf_and_conditions(samplesheet):
    for r in samplesheet:
        day = int(r["day_of_exposure"])
        assert int(r["dpf"]) == (19 if day == 0 else day + 18)
        assert r["condition"] == {"16": "baseline", "23": "control", "29": "treatment"}[r["temperature_C"]]
        assert r["raw_file"].endswith(".raw") and len(r["raw_sha1"]) == 40


def test_failed_qc_nothing_harmonized(study):
    assert study["qc_status"] == "failed"
    assert study["analysis_status"] == "not_started"
    (comp,) = study["comparisons"]
    assert comp["results_file"] is None
    assert comp["biological_replicates"] == 1
    assert comp["sample_size"] == 12
    assert not list(STUDY_DIR.glob("*_standardized.tsv"))
    for flag in ("low_replication", "tank_level_replication", "batch_confounded"):
        assert flag in study["quality_flags"]
    assert "No biological replication" in study["limitations"]


def test_phenotype_is_exposure_only(study):
    (comp,) = study["comparisons"]
    assert comp["comparison_id"] == "temp29_vs_temp23_days3to13"
    assert comp["meta_analysis_primary"] is True
    assert comp["stressor_standardized"] == "temperature"
    assert comp["phenotype"] == "growth_under_stress"
    assert comp["resilience_classification"] == "exposure_only"
    assert comp["control_condition"].startswith("23 C, Silo 3")
    assert comp["treatment_condition"].startswith("29 C, Silo 9")


def test_processed_table_is_located_but_not_committed(study, provenance):
    assert study["data_availability"]["processed_data_available"] is True
    assert provenance["figshare"]["md5"] == "d728060b4cac390fdeb7148c2dd772d3"
    assert provenance["abacus_output"] == {"n_proteins": 8443, "n_adjnsaf_columns": 26}
    assert provenance["search_database_deposited_in_pride"] is False
    assert provenance["search_database_in_figshare"] is True
    assert not list(STUDY_DIR.glob("*.zip")) and not (STUDY_DIR / "Abacus_output.tsv").exists()


def test_samplesheet_checksum_matches_provenance(provenance):
    rel = "data/studies/WANAMAKER2020_HEAT_SEED/samplesheet.tsv"
    digest = hashlib.sha256((REPO_ROOT / rel).read_bytes()).hexdigest()
    assert provenance["outputs"][rel] == digest
