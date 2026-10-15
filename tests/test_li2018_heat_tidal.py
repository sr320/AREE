"""Curation tests for LI2018_HEAT_TIDAL.

Li et al. 2018 (PRJNA407831): 35 C heat in adult gill of common-garden F1
from intertidal (tolerant) and subtidal (sensitive) populations at two
sites, one pooled library per tank at 0, 6 and 24 h. These tests pin the
design parsed from the aliases (cross-checked against the deposited `temp`
and `isolate` attributes), the authors' exclusion, the prespecified
comparisons, and that the tolerance contrast is fitted on populations.
"""
from __future__ import annotations

import csv
import json

import pytest
import yaml

from common import DATA_DIR, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "LI2018_HEAT_TIDAL.yaml"
STUDY_DIR = DATA_DIR / "studies" / "LI2018_HEAT_TIDAL"
POPULATIONS = ("BYQJ", "BYQX", "LTJ", "LTX")


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


def test_design_is_four_populations_by_three_times_by_three_tanks(samplesheet):
    prov = json.loads((STUDY_DIR / "ena_provenance.json").read_text())
    assert prov["bioproject"] == "PRJNA407831"
    assert prov["n_runs"] == len(samplesheet) == 36
    assert set(prov["group_sizes"].values()) == {3} and len(prov["group_sizes"]) == 12
    assert {r["population"] for r in samplesheet} == set(POPULATIONS)


def test_alias_parse_agrees_with_deposited_attributes(samplesheet):
    """Hours come from the alias but must match BioSample `temp` (35_<h>);
    the site must match `isolate`."""
    for row in samplesheet:
        assert row["temp"] == f"35_{row['hours']}"
        assert row["isolate"] == row["site"]
        assert row["population"] == row["site"] + row["zone"]


def test_comparisons_and_primary(study):
    ids = [c["comparison_id"] for c in study["comparisons"]]
    assert ids == ["heat_35c_24h_vs_0h", "heat_35c_6h_vs_0h", "intertidal_vs_subtidal_0h"]
    primary = [c["comparison_id"] for c in study["comparisons"] if c.get("meta_analysis_primary")]
    assert primary == ["heat_35c_24h_vs_0h"]


def test_heat_is_stress_response_and_tidal_origin_is_resilience(study):
    by_id = {c["comparison_id"]: c for c in study["comparisons"]}
    for cid in ("heat_35c_24h_vs_0h", "heat_35c_6h_vs_0h"):
        assert by_id[cid]["resilience_classification"] == "stress_response"
        assert by_id[cid]["sample_size"] == 23
    tidal = by_id["intertidal_vs_subtidal_0h"]
    assert tidal["resilience_classification"] == "resilience"
    assert tidal["phenotype_direction"] == "increase"
    # Replicates are populations, not tanks.
    assert tidal["biological_replicates"] == 2
    assert "family_level_replication" in study["quality_flags"]


def test_the_authors_outlier_is_excluded_everywhere(study, samplesheet):
    lt_j_0h_2 = [r for r in samplesheet if r["population"] == "LTJ" and r["hours"] == "0" and r["tank"] == "2"]
    assert [r["sample_id"] for r in lt_j_0h_2] == ["LTJ02"]
    for comp in study["comparisons"]:
        assert [e["sample_id"] for e in comp["excluded_samples"]] == ["LTJ02"]


def test_fastq_manifest_covers_every_run_paired(samplesheet):
    manifest = _tsv(STUDY_DIR / "fastq_manifest.tsv")
    assert {r["run_accession"] for r in manifest} == {r["run_accession"] for r in samplesheet}
    assert len(manifest) == 72
    assert all(len(r["md5"]) == 32 for r in manifest)


def test_run_outputs_match_the_registered_comparisons(study):
    from common import REPO_ROOT

    expected = {
        "heat_35c_24h_vs_0h": (["BYQ_J_0", "BYQ_X_0", "LT_J_0", "LT_X_0"], ["BYQ_J_24", "BYQ_X_24", "LT_J_24", "LT_X_24"], None),
        "heat_35c_6h_vs_0h": (["BYQ_J_0", "BYQ_X_0", "LT_J_0", "LT_X_0"], ["BYQ_J_6", "BYQ_X_6", "LT_J_6", "LT_X_6"], None),
        "intertidal_vs_subtidal_0h": (["BYQ_X_0", "LT_X_0"], ["BYQ_J_0", "LT_J_0"], "population"),
    }
    for comp in study["comparisons"]:
        cid = comp["comparison_id"]
        params = json.loads((STUDY_DIR / f"LI2018_HEAT_TIDAL_{cid}_workflow_manifest.json").read_text())["parameters"]
        control, treatment, unit = expected[cid]
        assert params["comparison_id"] == cid
        assert params["control_levels"] == control and params["treatment_levels"] == treatment
        assert params["replicate_unit"] == unit
        assert params["exclude_samples"] == ["LTJ02"]
        rows = _tsv(REPO_ROOT / comp["results_file"])
        assert len(rows) > 20000


def test_genotype_check_shows_distinct_animals_within_groups():
    text = (STUDY_DIR / "replicate_genotype_check.txt").read_text()
    within = [line for line in text.splitlines() if line.strip().startswith("within")][0]
    observed, noise = float(within.split("observed")[1].split()[0]), float(within.split("noise")[1].split()[0])
    assert observed > 2 * noise
