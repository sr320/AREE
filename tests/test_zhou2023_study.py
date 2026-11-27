"""Curation tests for ZHOU2023_SALINITY_GILL.

The first real hard clam (Mercenaria mercenaria) study: gill of clams held at
30 psu (NS), long-term 10 psu (CL) or long-term 45 psu (CH), three libraries
per group, each a pool of three clams. Four acute groups are registered in the
sample sheet but not analysed. These tests pin the design decoded from the
sample aliases, the pooled-library replication, the prespecified primary
comparison (hyposalinity), the exposure-only classification, the NCBI
download provenance for the nine analysed runs, and that both committed runs
are the registered ones (full-length path, 3 vs 3 pools, reference as fetched).
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_ID = "ZHOU2023_SALINITY_GILL"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
CONFIG_PATH = REPO_ROOT / "config" / f"{STUDY_ID}.config"
PRIMARY = "chronic_hyposal_10_vs_30psu"
SECONDARY = "chronic_hypersal_45_vs_30psu"
LEVELS = {PRIMARY: ("NS", "CL"), SECONDARY: ("NS", "CH")}
ANALYSED = {"NS", "CL", "CH"}
SALINITY = {
    "NS": "30 psu", "CL": "10 psu", "AL1": "10 psu", "AL5": "10 psu",
    "CH": "45 psu", "AH1": "45 psu", "AH5": "45 psu",
}


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


def _manifest(comparison_id):
    path = STUDY_DIR / f"{STUDY_ID}_{comparison_id}_workflow_manifest.json"
    return json.loads(path.read_text())


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_design_is_decoded_from_the_sample_aliases(samplesheet, ena_provenance):
    assert ena_provenance["bioproject"] == "PRJNA975821"
    assert len(samplesheet) == ena_provenance["n_runs"] == ena_provenance["n_samples"] == 21
    assert {r["library_layout"] for r in samplesheet} == {"PAIRED"}
    for row in samplesheet:
        m = re.fullmatch(r"(NS|CL|CH|AL1|AL5|AH1|AH5)_([123])", row["sample_alias"])
        assert m, row["sample_alias"]
        assert row["condition"] == m.group(1)
        assert row["store_cond"].startswith(SALINITY[m.group(1)])
        assert row["temp"] == "20℃"


def test_three_distinct_pools_per_group(samplesheet):
    assert Counter(r["condition"] for r in samplesheet) == dict.fromkeys(SALINITY, 3)
    assert len({r["sample_accession"] for r in samplesheet}) == 21
    assert {r["sieving"] for r in samplesheet} == {"3 individuals pooled into one sample"}


def test_hyposalinity_is_primary_and_both_are_exposure_only(study):
    comps = {c["comparison_id"]: c for c in study["comparisons"]}
    assert set(comps) == {PRIMARY, SECONDARY}
    assert comps[PRIMARY]["meta_analysis_primary"] is True
    assert not comps[SECONDARY].get("meta_analysis_primary", False)
    for comp in comps.values():
        assert comp["tissue"] == "gill"
        assert comp["resilience_classification"] == "exposure_only"
        assert comp["contrast_type"] == "exposure_response"
        assert comp["phenotype"] == "salinity_tolerance"
        assert comp["stressor_standardized"] == "salinity"
        assert comp["sample_size"] == 6 and comp["biological_replicates"] == 3


def test_pooling_is_recorded_as_a_limitation(study):
    assert study["species"] == "Mercenaria mercenaria"
    assert "Pooled libraries" in study["limitations"]


def test_both_comparisons_have_results(study):
    for comp in study["comparisons"]:
        assert (REPO_ROOT / comp["results_file"]).exists()
    assert study["analysis_status"] == "harmonized"
    assert study["qc_status"] == "passed_with_warnings"


@pytest.mark.parametrize("comparison_id", [PRIMARY, SECONDARY])
def test_the_committed_runs_are_the_ones_registered(comparison_id):
    """Full-length path, no covariates, 3 vs 3 pools, reference as fetched."""
    manifest = _manifest(comparison_id)
    control, treatment = LEVELS[comparison_id]
    assert manifest["parameters"]["tag_seq"] is False
    assert manifest["parameters"]["control_levels"] == [control]
    assert manifest["parameters"]["treatment_levels"] == [treatment]
    assert manifest["deseq2_design"]["design"] == "~condition"
    assert manifest["deseq2_design"]["length_model"].startswith("tximport")
    assert manifest["deseq2_design"]["n_libraries_or_units"] == 6
    assert manifest["deseq2_design"]["residual_df"] == 4
    assert all(f.get("matches_bundle_provenance", True) for f in manifest["reference"]["files"])
    assert manifest["warnings"] == []


def test_every_analysed_fastq_was_read_count_verified_from_ncbi(samplesheet):
    """The 9 NS, CL and CH runs (18 files) came from NCBI; each has a SHA-256."""
    with open(STUDY_DIR / "fastq_ncbi_provenance.tsv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    runs = {r["run_accession"] for r in samplesheet if r["condition"] in ANALYSED}
    assert len(runs) == 9
    assert len(rows) == 18
    assert {(r["run_accession"], r["mate"]) for r in rows} == {
        (run, mate) for run in runs for mate in ("1", "2")
    }
    for row in rows:
        assert row["source"] == "NCBI SRA"
        assert row["read_count_observed"] == row["read_count_expected"]
        assert re.fullmatch(r"[0-9a-f]{64}", row["sha256"])


def test_replicate_dispersion_qc_passes():
    with open(STUDY_DIR / "replicate_dispersion_qc.tsv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["condition"] for r in rows} == ANALYSED
    assert all(r["passes"] == "True" for r in rows)


def test_undeposited_exposure_details_stay_missing(study):
    for comp in study["comparisons"]:
        assert comp["exposure_duration"] is None and comp["exposure_timing"] is None
    assert {"low_replication", "ambiguous_phenotype_definition"} <= set(study["quality_flags"])


def test_config_runs_the_full_length_path_with_hyposalinity_as_default():
    text = CONFIG_PATH.read_text(encoding="utf-8")
    assert not re.search(r"^\s*tag_seq\s*=\s*true\b", text, re.M)
    assert not re.search(r"^\s*single_end\s*=\s*true\b", text, re.M)
    assert re.search(r"^\s*control_level\s*=\s*'NS'", text, re.M)
    assert re.search(r"^\s*treatment_level\s*=\s*'CL'", text, re.M)
    block = re.search(r"rnaseq \{(.*?)\}", text, re.S).group(1)
    assert re.search(r"comparison_id\s*=\s*null", block)
