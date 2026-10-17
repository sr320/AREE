"""Curation tests for PERUZZA2023_HEATWAVE.

The second real Manila clam study: 30 vs 20 C for 30 days, 20 clams per arm
from 4 tanks each, 3' Tag-Seq of digestive gland. These tests pin the design
read from ENA and the decisions that are easy to erode later: the 16S amplicon
runs that share the BioProject are excluded, the one prespecified comparison
and its exposure-only classification, the tank caveat, sex as a covariate, and
that the run took the workflow's 3' Tag-seq path (tag_seq = true).
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_ID = "PERUZZA2023_HEATWAVE"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
CONFIG_PATH = REPO_ROOT / "config" / f"{STUDY_ID}.config"
CID = "heatwave_30c_vs_control_20c"


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
    return json.loads((STUDY_DIR / f"{STUDY_ID}_{CID}_workflow_manifest.json").read_text())


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_samplesheet_keeps_only_the_rnaseq_runs(samplesheet, ena_provenance):
    """The BioProject also holds 39 paired-end 16S amplicon runs from the same RNA."""
    assert ena_provenance["bioproject"] == "PRJNA906528"
    assert ena_provenance["library_strategy_filter"] == {
        "kept": ["RNA-Seq"], "runs_by_strategy_before_filter": {"AMPLICON": 39, "RNA-Seq": 40}}
    assert len(samplesheet) == ena_provenance["n_runs"] == ena_provenance["n_samples"] == 40
    assert {r["library_layout"] for r in samplesheet} == {"SINGLE"}


def test_twenty_distinct_clams_per_temperature_with_unbalanced_sex(samplesheet):
    assert Counter(r["condition"] for r in samplesheet) == {"20": 20, "30": 20}
    assert len({r["sample_accession"] for r in samplesheet}) == 40
    assert len({r["sample_alias"] for r in samplesheet}) == 40
    sex = Counter((r["condition"], r["sex"]) for r in samplesheet)
    assert sex == {("20", "male"): 15, ("20", "female"): 5, ("30", "male"): 10, ("30", "female"): 10}


def test_one_primary_comparison_classified_exposure_only(study):
    """A non-lethal heatwave with no survival difference: exposure, not tolerance."""
    (comp,) = study["comparisons"]
    assert comp["comparison_id"] == CID and comp["meta_analysis_primary"] is True
    assert comp["resilience_classification"] == "exposure_only"
    assert comp["phenotype"] == "thermal_tolerance"
    assert comp["tissue"] == "digestive_gland" and comp["stressor_standardized"] == "temperature"
    assert comp["biological_replicates"] == 4


def test_tank_replication_tag_seq_and_status_are_declared(study):
    assert {"tank_level_replication", "unbalanced_design"} <= set(study["quality_flags"])
    assert "3' Tag" in study["platform"]
    assert study["analysis_status"] == "harmonized"
    assert study["qc_status"] == "passed_with_warnings"


def test_config_runs_the_tag_seq_path_with_sex_as_covariate(samplesheet):
    text = CONFIG_PATH.read_text(encoding="utf-8")
    assert re.search(r"^\s*tag_seq\s*=\s*true\b", text, re.M)
    assert re.search(r"^\s*covariates\s*=\s*'sex'", text, re.M)
    levels = set(re.findall(r"(?:control_level|treatment_level)\s*=\s*'([^']+)'", text))
    assert levels == {"20", "30"} <= {r["condition"] for r in samplesheet}
    block = re.search(r"rnaseq \{(.*?)\}", text, re.S).group(1)
    assert re.search(r"comparison_id\s*=\s*null", block)


def test_the_committed_run_is_the_one_registered(study, manifest):
    """tag_seq, no length offsets, sex fitted, 20 vs 20 clams, reference as fetched."""
    assert (REPO_ROOT / study["comparisons"][0]["results_file"]).exists()
    assert manifest["parameters"]["tag_seq"] is True
    assert manifest["deseq2_design"]["design"] == "~sex + condition"
    assert manifest["deseq2_design"]["length_model"].startswith("none")
    assert manifest["deseq2_design"]["n_libraries_or_units"] == 40
    assert manifest["deseq2_design"]["residual_df"] == 37
    assert all(f.get("matches_bundle_provenance", True) for f in manifest["reference"]["files"])
    assert manifest["warnings"] == []


def test_replicate_dispersion_qc_passes():
    with open(STUDY_DIR / "replicate_dispersion_qc.tsv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["condition"] for r in rows} == {"20", "30"}
    assert all(r["passes"] == "True" for r in rows)
