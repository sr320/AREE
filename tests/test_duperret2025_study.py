"""Curation tests for DUPERRET2025_POMS_STARVE.

Real metabolomics study (MetaboLights MTBLS13313; Duperret et al. 2025,
iScience) registered from deposited ISA-Tab metadata and the paper before
any analysis. These tests pin the design (2 families x 2 diets x 4 times x 3
tank pools), the finding that the deposited MAFs hold no abundance values,
the diet/MS-batch confounding, and the curation decisions that are easy to
erode later: tank-level replication, the single 0 hpc comparison the paper
reports, and the linked RNA-seq projects kept out of `accessions`.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "DUPERRET2025_POMS_STARVE.yaml"
STUDY_DIR = DATA_DIR / "studies" / "DUPERRET2025_POMS_STARVE"
PLATFORMS = ("GM", "NUCs", "PUFAs", "TCA")


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def samplesheet():
    with open(STUDY_DIR / "samplesheet.tsv") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def provenance():
    return json.loads((STUDY_DIR / "metabolights_provenance.json").read_text())


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_design_is_families_by_diets_by_times_by_three_tanks(samplesheet, provenance):
    assert provenance["accession"] == "MTBLS13313"
    assert provenance["n_samples"] == len(samplesheet) == 48
    cells = Counter((r["family"], r["condition"], r["timepoint_hpc"]) for r in samplesheet)
    assert set(cells) == {(f, d, t) for f in ("F11N", "F14R") for d in ("fed", "starved")
                          for t in ("0", "3", "6", "12")}
    assert set(cells.values()) == {3}
    for r in samplesheet:
        assert r["diet_code"] == {"fed": "3", "starved": "5"}[r["condition"]]
        assert r["sample_id"] == f"{r['family']}-{r['diet_code']}-{r['timepoint_hpc']}-{r['replicate']}"
        assert r["sample_unit"].startswith("pool of 10")


def test_every_sample_has_one_raw_file_per_platform(samplesheet):
    for platform in PLATFORMS:
        files = [r[f"raw_file_{platform}"] for r in samplesheet]
        assert all(f.startswith("FILES/RAW_FILES/") for f in files)
        assert len(set(files)) == 48


def test_deposited_mafs_hold_no_abundances(provenance, study):
    checks = provenance["maf_checks"]
    assert {p: c["n_features"] for p, c in checks.items()} == {"GM": 111, "NUCs": 21, "PUFAs": 183, "TCA": 13}
    for c in checks.values():
        assert c["n_abundance_values"] == 0
        assert c["n_sample_columns_in_samplesheet"] == 48
    assert study["data_availability"]["processed_data_available"] is False
    assert study["comparisons"][0]["results_file"] is None


def test_diet_is_confounded_with_ms_batch_and_says_so(provenance, study):
    batches = provenance["ms_batches_by_diet"]
    assert set(batches) == set(PLATFORMS)
    assert not any(b["diet_batches_overlap"] for b in batches.values())
    assert "confounded with diet" in study["limitations"]


def test_one_tank_replicated_0hpc_comparison(study):
    (comp,) = study["comparisons"]
    assert comp["comparison_id"] == "starved_vs_fed_0hpc"
    assert comp["meta_analysis_primary"] is True
    assert comp["stressor_standardized"] == "nutritional_limitation"
    assert comp["phenotype"] == "disease_resistance"
    assert comp["phenotype_direction"] == "increase"
    assert comp["resilience_classification"] == "resilience"
    assert comp["sample_size"] == 12 and comp["biological_replicates"] == 3
    assert comp["control_condition"].startswith("Fed ad libitum")
    for flag in ("tank_level_replication", "low_replication"):
        assert flag in study["quality_flags"]


def test_linked_sequencing_projects_are_context_not_accessions(study):
    accs = " ".join(str(v) for v in study["accessions"].values() if v)
    for acc in ("PRJEB86618", "PRJEB86525", "PRJEB86647"):
        assert acc not in accs
        assert acc in study["citation"]
    assert study["accessions"]["doi"] == "10.1016/j.isci.2025.114333"
    assert study["accessions"]["other"] == "MTBLS13313"


def test_samplesheet_checksum_matches_provenance(provenance):
    rel = "data/studies/DUPERRET2025_POMS_STARVE/samplesheet.tsv"
    digest = hashlib.sha256((REPO_ROOT / rel).read_bytes()).hexdigest()
    assert provenance["outputs"][rel] == digest
