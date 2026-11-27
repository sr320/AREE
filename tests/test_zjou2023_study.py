"""Curation tests for ZJOU2023_HEAT.

The first real Mediterranean mussel (Mytilus galloprovincialis) study: mussels
from Zhoushan held at 30 C vs 18 C, gill and mantle, three mussels per tissue
and temperature. These tests pin the design decoded from the sample aliases,
the prespecified primary comparison (gill), the exposure-only classification,
the NCBI download provenance, and that both committed runs are the registered
ones (full-length path, 3 vs 3 mussels, one tissue each, reference as fetched).
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_ID = "ZJOU2023_HEAT"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
CONFIG_PATH = REPO_ROOT / "config" / f"{STUDY_ID}.config"
PRIMARY = "heat_30c_vs_18c_gill"
SECONDARY = "heat_30c_vs_18c_mantle"
LEVELS = {PRIMARY: ("G_C", "G_H"), SECONDARY: ("M_C", "M_H")}
TISSUE = {PRIMARY: "gill", SECONDARY: "mantle"}
TEMP = {"C": "18℃", "H": "30℃"}


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
    assert ena_provenance["bioproject"] == "PRJNA954971"
    assert len(samplesheet) == ena_provenance["n_runs"] == ena_provenance["n_samples"] == 12
    assert {r["library_layout"] for r in samplesheet} == {"PAIRED"}
    for row in samplesheet:
        m = re.fullmatch(r"([GM])([CH])([123])", row["sample_alias"])
        assert m, row["sample_alias"]
        tissue_code, temp_code, _ = m.groups()
        assert row["condition"] == f"{tissue_code}_{temp_code}"
        assert row["temp"] == TEMP[temp_code]


def test_three_distinct_mussels_per_tissue_and_temperature(samplesheet):
    assert Counter(r["condition"] for r in samplesheet) == {
        "G_C": 3, "G_H": 3, "M_C": 3, "M_H": 3,
    }
    assert len({r["sample_accession"] for r in samplesheet}) == 12


def test_gill_is_the_only_pooled_comparison_and_both_are_exposure_only(study):
    comps = {c["comparison_id"]: c for c in study["comparisons"]}
    assert set(comps) == {PRIMARY, SECONDARY}
    assert comps[PRIMARY]["meta_analysis_primary"] is True
    assert not comps[SECONDARY].get("meta_analysis_primary", False)
    for cid, comp in comps.items():
        assert comp["tissue"] == TISSUE[cid]
        assert comp["resilience_classification"] == "exposure_only"
        assert comp["contrast_type"] == "exposure_response"
        assert comp["phenotype"] == "thermal_tolerance"
        assert comp["stressor_standardized"] == "temperature"
        assert comp["sample_size"] == 6 and comp["biological_replicates"] == 3


def test_species_label_is_flagged_as_unverified(study):
    assert study["species"] == "Mytilus galloprovincialis"
    assert study["species_identity"]["basis"] == "authors_statement"
    assert "species_assignment_unverified" in study["quality_flags"]


def test_both_comparisons_have_results(study):
    for comp in study["comparisons"]:
        assert (REPO_ROOT / comp["results_file"]).exists()
    assert study["analysis_status"] == "harmonized"
    assert study["qc_status"] == "passed_with_warnings"


@pytest.mark.parametrize("comparison_id", [PRIMARY, SECONDARY])
def test_the_committed_runs_are_the_ones_registered(comparison_id):
    """Full-length path, no covariates, 3 vs 3 mussels of one tissue, reference as fetched."""
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


def test_every_fastq_was_read_count_verified_from_ncbi(samplesheet):
    """All 24 files came from NCBI fastq-dump; each has a SHA-256."""
    with open(STUDY_DIR / "fastq_ncbi_provenance.tsv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    runs = {r["run_accession"] for r in samplesheet}
    assert len(rows) == 24
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
    assert {r["condition"] for r in rows} == {"G_C", "G_H", "M_C", "M_H"}
    assert all(r["passes"] == "True" for r in rows)


def test_undeposited_exposure_details_stay_missing(study):
    for comp in study["comparisons"]:
        assert comp["exposure_duration"] is None and comp["exposure_timing"] is None
    assert {"low_replication", "ambiguous_phenotype_definition"} <= set(study["quality_flags"])
    assert study["accessions"]["doi"] is None


def test_config_runs_the_full_length_path_with_gill_as_default():
    text = CONFIG_PATH.read_text(encoding="utf-8")
    assert not re.search(r"^\s*tag_seq\s*=\s*true\b", text, re.M)
    assert not re.search(r"^\s*single_end\s*=\s*true\b", text, re.M)
    assert re.search(r"^\s*control_level\s*=\s*'G_C'", text, re.M)
    assert re.search(r"^\s*treatment_level\s*=\s*'G_H'", text, re.M)
    block = re.search(r"rnaseq \{(.*?)\}", text, re.S).group(1)
    assert re.search(r"comparison_id\s*=\s*null", block)
