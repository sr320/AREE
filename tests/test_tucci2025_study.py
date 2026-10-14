"""Curation tests for TUCCI2025_HEAT_PRIMING.

The first real clam study: a 2x2 Manila clam priming x heatwave factorial,
10 clams per group, 3' Tag-Seq of digestive gland. These tests pin the design
read from ENA and the decisions that are easy to erode later: the 16S amplicon
runs that share the BioProject are excluded, the two prespecified comparisons
and their classifications, the tank and survivor caveats, and that the study is
configuration runs the workflow's 3' Tag-seq path (tag_seq = true).
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_ID = "TUCCI2025_HEAT_PRIMING"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
CONFIG_PATH = REPO_ROOT / "config" / f"{STUDY_ID}.config"
LEVELS = {"Naïve_No_Heatwave", "Naïve_Heatwave", "Primed_No_Heatwave", "Primed_Heatwave"}


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


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_samplesheet_keeps_only_the_rnaseq_runs(samplesheet, ena_provenance):
    """The BioProject also holds 40 paired-end 16S amplicon runs from the same RNA.
    Counting them would double every group."""
    assert ena_provenance["bioproject"] == "PRJNA1132140"
    assert ena_provenance["library_strategy_filter"] == {
        "kept": ["RNA-Seq"], "runs_by_strategy_before_filter": {"AMPLICON": 40, "RNA-Seq": 40}}
    assert ena_provenance["library_layouts"] == ["SINGLE"]
    assert len(samplesheet) == ena_provenance["n_runs"] == ena_provenance["n_samples"] == 40
    assert {r["library_layout"] for r in samplesheet} == {"SINGLE"}


def test_ten_distinct_clams_per_group(samplesheet):
    assert Counter(r["condition"] for r in samplesheet) == dict.fromkeys(LEVELS, 10)
    assert len({r["clam_id"] for r in samplesheet}) == 40
    assert len({r["sample_accession"] for r in samplesheet}) == 40


def test_two_prespecified_comparisons_with_one_primary(study):
    comparisons = {c["comparison_id"]: c for c in study["comparisons"]}
    assert set(comparisons) == {"naive_hw_vs_naive_control", "primed_hw_vs_naive_hw"}
    primary = [cid for cid, c in comparisons.items() if c.get("meta_analysis_primary")]
    assert primary == ["naive_hw_vs_naive_control"]


def test_classifications_follow_what_was_measured(study):
    """Heatwave vs control is an exposure among survivors. Primed vs naive under
    the same heatwave compares groups with a measured survival difference."""
    comparisons = {c["comparison_id"]: c for c in study["comparisons"]}
    assert comparisons["naive_hw_vs_naive_control"]["resilience_classification"] == "exposure_only"
    assert comparisons["primed_hw_vs_naive_hw"]["resilience_classification"] == "resilience"
    for c in comparisons.values():
        assert c["tissue"] == "digestive_gland"
        assert c["stressor_standardized"] == "temperature"
        assert c["phenotype"] == "thermal_tolerance"
        assert "survivors only" in c["exposure_timing"]


def test_tank_replication_and_tag_seq_are_declared(study):
    assert "tank_level_replication" in study["quality_flags"]
    assert all(c["biological_replicates"] == 2 for c in study["comparisons"])
    assert "3' Tag" in study["platform"]
    assert study["analysis_status"] == "harmonized"
    assert study["qc_status"] == "passed_with_warnings"


def test_config_levels_exist_in_the_samplesheet(samplesheet):
    text = CONFIG_PATH.read_text(encoding="utf-8")
    levels = set(re.findall(r"(?:control_level|treatment_level)\s*=\s*'([^']+)'", text))
    assert levels == {"Naïve_No_Heatwave", "Naïve_Heatwave"}
    assert levels <= {r["condition"] for r in samplesheet}
    assert "single_end = true" in text
    assert "GCF_026571515.1" in text


def test_config_runs_the_tag_seq_path():
    """QuantSeq 3' Tag-Seq: the standard path's length offsets would create false
    fold changes from isoform switching. Never run this study without tag_seq."""
    text = CONFIG_PATH.read_text(encoding="utf-8")
    assert re.search(r"^\s*tag_seq\s*=\s*true\b", text, re.M)


def test_command_line_comparison_id_is_not_shadowed():
    """The workflow reads rnaseq.comparison_id before comparison_id; set inside the
    rnaseq block it would silently name the second comparison after the first."""
    text = CONFIG_PATH.read_text(encoding="utf-8")
    block = re.search(r"rnaseq \{(.*?)\}", text, re.S).group(1)
    assert re.search(r"comparison_id\s*=\s*null", block)
    assert re.search(r"^\s*comparison_id\s*=\s*'naive_hw_vs_naive_control'", text, re.M)


@pytest.mark.parametrize("cid", ["naive_hw_vs_naive_control", "primed_hw_vs_naive_hw"])
def test_each_comparison_was_run_on_the_tag_seq_path(study, cid):
    """The committed manifest must show the run that produced the results file:
    tag_seq, no length offsets, 10 vs 10 clams, reference files as fetched."""
    comp = {c["comparison_id"]: c for c in study["comparisons"]}[cid]
    assert (REPO_ROOT / comp["results_file"]).exists()
    manifest = json.loads((STUDY_DIR / f"{STUDY_ID}_{cid}_workflow_manifest.json").read_text())
    assert manifest["parameters"]["tag_seq"] is True
    assert manifest["deseq2_design"]["length_model"].startswith("none")
    assert manifest["deseq2_design"]["n_libraries_or_units"] == 20
    assert all(f.get("matches_bundle_provenance", True) for f in manifest["reference"]["files"])
    assert manifest["warnings"] == []


def test_replicate_dispersion_qc_passes():
    with open(STUDY_DIR / "replicate_dispersion_qc.tsv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["condition"] for r in rows} == LEVELS
    assert all(r["passes"] == "True" for r in rows)
