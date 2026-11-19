"""Curation tests for NORTHSOUTH2024_TEMP_GILL.

The third real Manila clam study: gill of two wild populations (North,
Golden Pebble Beach, Dalian; South, Guangxi) held at -1, 12 or 28 C, three
clams per cell. These tests pin the design decoded from the sample aliases,
the prespecified primary (heat) comparison with population as a covariate,
the cold comparison kept out of pooling, the exposure-only classification,
and that the committed heat run is the registered one (full-length path,
population fitted, 6 vs 6 clams, reference as fetched).
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_ID = "NORTHSOUTH2024_TEMP_GILL"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
CONFIG_PATH = REPO_ROOT / "config" / f"{STUDY_ID}.config"
PRIMARY = "heat_28c_vs_12c"
COLD = "cold_minus1c_vs_12c"


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def samplesheet():
    with open(STUDY_DIR / "samplesheet.tsv", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def ena_provenance():
    return json.loads((STUDY_DIR / "ena_provenance.json").read_text())


@pytest.fixture(scope="module")
def manifest():
    return json.loads((STUDY_DIR / f"{STUDY_ID}_{PRIMARY}_workflow_manifest.json").read_text())


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_design_is_decoded_from_the_sample_aliases(samplesheet, ena_provenance):
    """G/J = South/North, 1/2/3 = -1/12/28 C; the alias must agree with the site."""
    assert ena_provenance["bioproject"] == "PRJNA738278"
    assert ena_provenance["run_labels"]["field"] == "sample_alias"
    assert len(samplesheet) == ena_provenance["n_runs"] == ena_provenance["n_samples"] == 18
    assert {r["library_layout"] for r in samplesheet} == {"PAIRED"}
    site = {"G": "Guangxi province", "J": "Golden Pebble Beach"}
    for row in samplesheet:
        assert row["sample_alias"][0] == row["population"]
        assert row["isolation_source"] == site[row["population"]]
        assert row["condition"] == row["temp"] == row["sample_alias"][1]


def test_three_distinct_clams_per_population_and_temperature(samplesheet):
    cells = Counter((r["population"], r["condition"]) for r in samplesheet)
    assert cells == {(p, t): 3 for p in "GJ" for t in "123"}
    assert len({r["sample_accession"] for r in samplesheet}) == 18


def test_heat_is_the_only_pooled_comparison_and_both_are_exposure_only(study):
    comps = {c["comparison_id"]: c for c in study["comparisons"]}
    assert set(comps) == {PRIMARY, COLD}
    assert comps[PRIMARY]["meta_analysis_primary"] is True
    assert not comps[COLD].get("meta_analysis_primary", False)
    for comp in comps.values():
        assert comp["resilience_classification"] == "exposure_only"
        assert comp["contrast_type"] == "exposure_response"
        assert comp["phenotype"] == "thermal_tolerance"
        assert comp["tissue"] == "gill" and comp["stressor_standardized"] == "temperature"
        assert comp["sample_size"] == 12 and comp["biological_replicates"] == 6


def test_only_the_heat_comparison_has_results(study):
    comps = {c["comparison_id"]: c for c in study["comparisons"]}
    assert (REPO_ROOT / comps[PRIMARY]["results_file"]).exists()
    assert comps[COLD]["results_file"] is None
    assert study["analysis_status"] == "harmonized"
    assert study["qc_status"] == "passed_with_warnings"


def test_the_committed_run_is_the_one_registered(manifest):
    """Full-length path, population fitted, 6 vs 6 clams, reference as fetched."""
    assert manifest["parameters"]["tag_seq"] is False
    assert manifest["parameters"]["control_levels"] == ["2"]
    assert manifest["parameters"]["treatment_levels"] == ["3"]
    assert manifest["deseq2_design"]["design"] == "~population + condition"
    assert manifest["deseq2_design"]["length_model"].startswith("tximport")
    assert manifest["deseq2_design"]["n_libraries_or_units"] == 12
    assert manifest["deseq2_design"]["residual_df"] == 9
    assert all(f.get("matches_bundle_provenance", True) for f in manifest["reference"]["files"])
    assert manifest["warnings"] == []


def test_replicate_dispersion_qc_passes():
    with open(STUDY_DIR / "replicate_dispersion_qc.tsv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["condition"] for r in rows} == {"2", "3"}
    assert all(r["passes"] == "True" for r in rows)


def test_undeposited_exposure_details_stay_missing(study):
    for comp in study["comparisons"]:
        assert comp["exposure_duration"] is None and comp["exposure_timing"] is None
    assert {"low_replication", "ambiguous_phenotype_definition"} <= set(study["quality_flags"])
    assert study["accessions"]["doi"] is None


def test_config_fits_population_and_runs_the_full_length_path(samplesheet):
    text = CONFIG_PATH.read_text(encoding="utf-8")
    assert not re.search(r"^\s*tag_seq\s*=\s*true\b", text, re.M)
    assert not re.search(r"^\s*single_end\s*=\s*true\b", text, re.M)
    assert re.search(r"^\s*covariates\s*=\s*'population'", text, re.M)
    assert re.search(r"^\s*control_level\s*=\s*'2'", text, re.M)
    assert re.search(r"^\s*treatment_level\s*=\s*'3'", text, re.M)
    assert {"2", "3"} <= {r["condition"] for r in samplesheet}
    block = re.search(r"rnaseq \{(.*?)\}", text, re.S).group(1)
    assert re.search(r"comparison_id\s*=\s*null", block)
