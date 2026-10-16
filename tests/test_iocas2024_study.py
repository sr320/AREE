"""Curation tests for IOCAS2024_OA_TRANSGEN.

Transgenerational ocean acidification in larvae (PRJNA1196326). It FAILED QC
before harmonization: each arm is one larval pool sequenced as three
libraries, and the arms descend from different parents (effective n=1 per
arm). These tests pin the design read from ENA and the decisions that are
easy to erode later: exposure_only classification, no borrowed DOI, and above
all that nothing from this study is harmonized or pooled.
"""
from __future__ import annotations

import csv
import json

import pytest
import yaml

from common import DATA_DIR, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "IOCAS2024_OA_TRANSGEN.yaml"
STUDY_DIR = DATA_DIR / "studies" / "IOCAS2024_OA_TRANSGEN"
COMPARISON = "parental_ph7_4_vs_ph8_1_larvae"


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


def test_samplesheet_matches_live_ena_design_snapshot(samplesheet, ena_provenance):
    assert ena_provenance["bioproject"] == "PRJNA1196326"
    assert ena_provenance["n_runs"] == ena_provenance["n_samples"] == len(samplesheet) == 6
    assert ena_provenance["group_sizes"] == {"74": 3, "81": 3}
    assert len({r["sample_alias"] for r in samplesheet}) == 6


def test_parental_ph_arms_agree_with_biosample_isolate(samplesheet):
    """The alias prefix (S74/S81) is the only machine-readable pH; check it
    against the free-text isolate so a relabelled deposit is caught."""
    for row in samplesheet:
        assert row["sample_alias"].startswith(f"S{row['condition']}_")
        expected = "acidification group" if row["condition"] == "74" else "wild type group"
        assert row["isolate"].startswith(expected)
        assert row["tissue"] == "larvae"


def test_no_measured_phenotype_means_exposure_only(study):
    (comp,) = study["comparisons"]
    assert comp["comparison_id"] == COMPARISON
    assert comp["stressor_standardized"] == "ocean_acidification"
    assert comp["resilience_classification"] == "exposure_only"
    assert comp["phenotype_direction"] == "not_applicable"
    assert comp["life_stage"] == "larva"
    assert comp["tissue"] == "whole_larva"
    # Parental exposure length and larval age are not deposited: keep them unknown.
    assert comp["exposure_duration"] is None
    assert "ambiguous_phenotype_definition" in study["quality_flags"]


def test_failed_qc_keeps_the_study_out_of_every_pool(study):
    """One pool per arm, one parental cohort per arm: nothing may be pooled."""
    assert study["qc_status"] == "failed"
    assert "low_replication" in study["quality_flags"]
    assert not any(c.get("meta_analysis_primary") for c in study["comparisons"])
    assert all(c["results_file"] is None for c in study["comparisons"])


def test_archived_dispersion_qc_records_the_failure():
    with open(STUDY_DIR / "replicate_dispersion_qc.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["condition"] for r in rows} == {"74", "81"}
    assert all(float(r["median_dispersion"]) < 0.005 for r in rows)
    assert all(r["passes"] == "False" for r in rows)


def test_run_is_archived_but_not_harmonized():
    manifest = STUDY_DIR / f"IOCAS2024_OA_TRANSGEN_{COMPARISON}_workflow_manifest.json"
    assert json.loads(manifest.read_text())["comparison_id"] == COMPARISON
    assert (STUDY_DIR / f"IOCAS2024_OA_TRANSGEN_{COMPARISON}_rnaseq_report.html").exists()
    assert not (STUDY_DIR / f"IOCAS2024_OA_TRANSGEN_{COMPARISON}_dge_standardized.tsv").exists()


def test_no_publication_is_attached(study):
    assert study["accessions"]["doi"] is None


def test_fastq_manifest_covers_every_run_as_paired_end(samplesheet):
    with open(STUDY_DIR / "fastq_manifest.tsv") as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    runs_in_sheet = {r["run_accession"] for r in samplesheet}
    assert {r["run_accession"] for r in manifest} == runs_in_sheet
    assert len(manifest) == 2 * len(runs_in_sheet)
    assert all(len(r["md5"]) == 32 for r in manifest)
