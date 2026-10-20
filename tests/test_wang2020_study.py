"""Curation tests for WANG2020_TIDAL_HEAT_WGBS.

WGBS of Laoting intertidal and subtidal F1 oysters before and after 6 h at
35 C: the methylation layer on LI2018_HEAT_TIDAL's Laoting animals. All 12
libraries share one BioSample, so the design is parsed from library_name,
which carries the paper's group codes. These tests pin that parse, the four
contrast types, the primary chosen among the two heat responses, and the
decisions that are easy to erode later: replication counted in populations
for the tidal-origin contrasts, and the study recorded as not independent of
LI2018.
"""
from __future__ import annotations

import csv
import json

import pytest
import yaml

from common import DATA_DIR, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "WANG2020_TIDAL_HEAT_WGBS.yaml"
STUDY_DIR = DATA_DIR / "studies" / "WANG2020_TIDAL_HEAT_WGBS"

# Run -> paper group, from ENA library_name (InC_1 ... SubH_3).
EXPECTED_GROUPS = {
    "SRR10065364": "SubH_1",
    "SRR10065365": "SubC_3",
    "SRR10065366": "SubC_2",
    "SRR10065367": "SubC_1",
    "SRR10065368": "InH_3",
    "SRR10065369": "InH_2",
    "SRR10065370": "InH_1",
    "SRR10065371": "InC_3",
    "SRR10065372": "SubH_3",
    "SRR10065373": "SubH_2",
    "SRR10065374": "InC_2",
    "SRR10065375": "InC_1",
}
# Wang et al. 2020 Table 1, mean raw reads per group.
PAPER_TABLE1_RAW_READS = {"In_C": 83021662, "In_H": 69956461, "Sub_C": 75452051, "Sub_H": 79864515}

EXPECTED_CONTRASTS = {
    "intertidal_vs_subtidal_0h": "baseline_phenotype_difference",
    "intertidal_vs_subtidal_6h_heat": "phenotype_difference_under_exposure",
    "heat_35c_6h_vs_0h_subtidal": "exposure_response",
    "heat_35c_6h_vs_0h_intertidal": "exposure_response",
}


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def comparisons(study):
    return {c["comparison_id"]: c for c in study["comparisons"]}


@pytest.fixture(scope="module")
def samplesheet():
    with open(STUDY_DIR / "samplesheet.tsv") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def ena_provenance():
    return json.loads((STUDY_DIR / "ena_provenance.json").read_text())


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_samplesheet_matches_the_deposit(samplesheet, ena_provenance):
    assert ena_provenance["bioproject"] == "PRJNA562805"
    assert ena_provenance["library_strategies"] == ["Bisulfite-Seq"]
    assert ena_provenance["n_runs"] == len(samplesheet) == 12
    # One BioSample for every library: the reason run labels were needed.
    assert ena_provenance["n_samples"] == 1
    assert ena_provenance["group_sizes"] == {"In_C": 3, "In_H": 3, "Sub_C": 3, "Sub_H": 3}
    assert ena_provenance["run_labels"]["field"] == "library_name"
    assert {r["run_accession"]: r["sample_id"] for r in samplesheet} == EXPECTED_GROUPS


def test_group_parse_is_consistent_with_labels(samplesheet):
    for row in samplesheet:
        assert row["sample_id"] == f"{row['origin']}{row['temperature']}_{row['replicate']}"
        assert row["condition"] == f"{row['origin']}_{row['temperature']}"
        assert row["tissue"] == "gill"


def test_group_read_counts_reproduce_paper_table1(samplesheet):
    """Independent check that library_name groups are the paper's groups."""
    for group, paper_mean in PAPER_TABLE1_RAW_READS.items():
        counts = [int(r["read_count"]) for r in samplesheet if r["condition"] == group]
        assert len(counts) == 3
        assert abs(sum(counts) / 3 - paper_mean) <= 1


def test_methylkit_sheet_codes_subtidal_as_baseline_reference(samplesheet):
    with open(STUDY_DIR / "methylkit_samplesheet.csv") as fh:
        coded = {r["sample_id"]: r["treatment"] for r in csv.DictReader(fh)}
    expected = {
        r["run_accession"]: "1" if r["condition"] == "In_C" else "0"
        for r in samplesheet
        if r["condition"] in ("In_C", "Sub_C")
    }
    assert coded == expected


def test_contrast_types(comparisons):
    assert {cid: c["contrast_type"] for cid, c in comparisons.items()} == EXPECTED_CONTRASTS


def test_tidal_origin_contrasts_are_resilience_with_one_population_per_class(comparisons):
    for cid in ("intertidal_vs_subtidal_0h", "intertidal_vs_subtidal_6h_heat"):
        comp = comparisons[cid]
        assert comp["resilience_classification"] == "resilience"
        assert comp["phenotype"] == "thermal_tolerance"
        assert comp["phenotype_direction"] == "increase"
        # Intertidal is the treatment arm, as in LI2018.
        assert comp["treatment_condition"].startswith("intertidal")
        assert comp["control_condition"].startswith("subtidal")
        assert comp["sample_size"] == 6
        assert comp["biological_replicates"] == 1


def test_heat_responses_are_stress_responses(comparisons):
    for cid in ("heat_35c_6h_vs_0h_subtidal", "heat_35c_6h_vs_0h_intertidal"):
        comp = comparisons[cid]
        assert comp["resilience_classification"] == "stress_response"
        assert comp["phenotype_direction"] == "not_applicable"
        assert comp["exposure_duration"] == "6 h"
        assert comp["sample_size"] == 6
        assert comp["biological_replicates"] == 3


def test_one_primary_per_shared_contrast_type(comparisons):
    primaries = {cid for cid, c in comparisons.items() if c.get("meta_analysis_primary")}
    assert primaries == {"heat_35c_6h_vs_0h_subtidal"}


def test_registered_with_run_pending(study, comparisons):
    assert study["analysis_mode"] == "raw_reanalysis"
    assert study["analysis_status"] == "not_started"
    assert study["qc_status"] == "not_started"
    assert all(c["results_file"] is None for c in comparisons.values())


def test_not_independent_of_li2018_is_recorded(study):
    assert "PRJNA407831" in study["accessions"]["other"]
    assert "NOT INDEPENDENT OF LI2018_HEAT_TIDAL" in study["limitations"]
    for flag in ("family_level_replication", "low_replication"):
        assert flag in study["quality_flags"]


def test_fastq_manifest_covers_every_run_as_paired_end(samplesheet):
    with open(STUDY_DIR / "fastq_manifest.tsv") as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    runs_in_sheet = {r["run_accession"] for r in samplesheet}
    assert {r["run_accession"] for r in manifest} == runs_in_sheet
    assert len(manifest) == 2 * len(runs_in_sheet)
    assert all(len(r["md5"]) == 32 for r in manifest)
