"""Curation tests for CIBNOR2021_HEAT_RRSS.

AREE's first real methylation study, and its first real thermal-resilience
contrast, registered from deposited metadata before any reanalysis. The
project deposits all 12 libraries under one BioSample, so the design is parsed
from experiment_title. These tests pin that parse, the family structure behind
it, and the curation decisions that are easy to erode later: biological
replicates counted as families, no borrowed DOI, and SS as the methylKit
reference.
"""
from __future__ import annotations

import csv
import json

import pytest
import yaml

from common import DATA_DIR, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "CIBNOR2021_HEAT_RRSS.yaml"
STUDY_DIR = DATA_DIR / "studies" / "CIBNOR2021_HEAT_RRSS"
FAMILIES = {"RR": {"52", "59"}, "SS": {"05", "35"}}


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


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


def test_design_is_two_families_per_class_three_libraries_each(samplesheet, ena_provenance):
    assert ena_provenance["bioproject"] == "PRJNA690951"
    assert ena_provenance["library_strategies"] == ["Bisulfite-Seq"]
    assert ena_provenance["n_runs"] == len(samplesheet) == 12
    # One BioSample for every library: the reason run labels were needed.
    assert ena_provenance["n_samples"] == 1
    assert ena_provenance["group_sizes"] == {"RR": 6, "SS": 6}

    for cls, fams in FAMILIES.items():
        rows = [r for r in samplesheet if r["condition"] == cls]
        assert {r["family"] for r in rows} == fams
        for fam in fams:
            assert sorted(r["replicate"] for r in rows if r["family"] == fam) == ["1", "2", "3"]


def test_parsed_labels_are_unique_and_recorded(samplesheet, ena_provenance):
    labels = ena_provenance["run_labels"]["labels"]
    assert ena_provenance["run_labels"]["field"] == "experiment_title"
    assert len(set(labels.values())) == 12
    for row in samplesheet:
        assert labels[row["run_accession"]] == row["sample_id"]
        assert row["sample_id"] == f"{row['condition']}{row['family']}-R{row['replicate']}"


def test_methylkit_sheet_codes_ss_as_reference_by_run_accession(samplesheet):
    with open(STUDY_DIR / "methylkit_samplesheet.csv") as fh:
        coded = {r["sample_id"]: r["treatment"] for r in csv.DictReader(fh)}
    expected = {r["run_accession"]: "1" if r["condition"] == "RR" else "0" for r in samplesheet}
    assert coded == expected


def test_replication_is_counted_in_families_not_libraries(study):
    (comp,) = study["comparisons"]
    assert comp["sample_size"] == 12
    assert comp["biological_replicates"] == 2
    assert comp["meta_analysis_primary"] is True
    for flag in ("family_level_replication", "low_replication", "ambiguous_phenotype_definition"):
        assert flag in study["quality_flags"]


def test_tolerance_contrast_under_one_challenge_not_an_exposure(study):
    (comp,) = study["comparisons"]
    assert comp["stressor_standardized"] == "temperature"
    assert comp["phenotype"] == "thermal_tolerance"
    assert comp["resilience_classification"] == "resilience"
    assert comp["phenotype_direction"] == "increase"
    assert comp["tissue"] == "gill"
    assert "not an unexposed control" in comp["control_condition"]
    assert comp["results_file"] is None


def test_rnaseq_paper_is_context_not_the_source_doi(study):
    assert study["accessions"]["doi"] is None
    assert "10.1016/j.cbd.2023.101089" in study["citation"]
    assert "does not analyze these data" in study["citation"]


def test_fastq_manifest_covers_every_run_as_paired_end(samplesheet):
    with open(STUDY_DIR / "fastq_manifest.tsv") as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    runs_in_sheet = {r["run_accession"] for r in samplesheet}
    assert {r["run_accession"] for r in manifest} == runs_in_sheet
    assert len(manifest) == 2 * len(runs_in_sheet)
    assert all(len(r["md5"]) == 32 for r in manifest)

