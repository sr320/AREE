"""Curation tests for IOCAS2022_OA_ENERGY.

AREE's first real ocean-acidification study. It was reanalyzed at full depth
and FAILED QC: its replicates are technical, not biological (effective n=1 per
group). These tests pin the design read from ENA and the decisions that are
easy to erode later: exposure_only classification, no borrowed DOI, and above
all that nothing from this study is harmonized or pooled.
"""
from __future__ import annotations

import csv
import json

import pytest
import yaml

from common import DATA_DIR, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "IOCAS2022_OA_ENERGY.yaml"
STUDY_DIR = DATA_DIR / "studies" / "IOCAS2022_OA_ENERGY"
TIMEPOINTS = ("7d", "28d", "56d")


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
    assert ena_provenance["bioproject"] == "PRJNA826964"
    assert ena_provenance["n_runs"] == ena_provenance["n_samples"] == len(samplesheet) == 18
    assert ena_provenance["min_group_size"] == 3
    # One library per animal: this is what rejected PRJNA623063 lacked.
    assert len({r["sample_alias"] for r in samplesheet}) == 18


def test_each_comparison_is_against_its_time_matched_control(study, samplesheet):
    by_condition: dict[str, int] = {}
    for row in samplesheet:
        by_condition[row["condition"]] = by_condition.get(row["condition"], 0) + 1

    ids = [c["comparison_id"] for c in study["comparisons"]]
    assert ids == [f"oa_{t}_vs_control_{t}" for t in TIMEPOINTS]
    for t in TIMEPOINTS:
        assert by_condition[f"control_{t}"] == 3
        assert by_condition[f"ocean_acidification_{t}"] == 3


def test_no_measured_phenotype_means_exposure_only(study):
    for comp in study["comparisons"]:
        assert comp["stressor_standardized"] == "ocean_acidification"
        assert comp["resilience_classification"] == "exposure_only"
        assert comp["phenotype_direction"] == "not_applicable"
        assert comp["tissue"] == "digestive_gland"
        assert comp["sample_size"] == 6
    assert "ambiguous_phenotype_definition" in study["quality_flags"]


def test_failed_qc_keeps_the_study_out_of_every_pool(study):
    """Pseudo-replicated (effective n=1): nothing may be harmonized or pooled."""
    assert study["qc_status"] == "failed"
    assert "low_replication" in study["quality_flags"]
    assert not any(c.get("meta_analysis_primary") for c in study["comparisons"])
    assert all(c["results_file"] is None for c in study["comparisons"])


def test_archived_dispersion_qc_records_the_failure():
    """Every group's replicates vary like technical replicates (~Poisson).

    Real M. gigas biological replicates sit at ~0.03-0.05 (CALLA2026_OSHV);
    this study's groups are one to two orders of magnitude lower.
    """
    with open(STUDY_DIR / "replicate_dispersion_qc.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 6
    assert all(float(r["median_dispersion"]) < 0.005 for r in rows)
    assert all(r["passes"] == "False" for r in rows)


def test_every_run_is_archived_but_none_is_harmonized():
    for t in TIMEPOINTS:
        cid = f"oa_{t}_vs_control_{t}"
        manifest = STUDY_DIR / f"IOCAS2022_OA_ENERGY_{cid}_workflow_manifest.json"
        assert json.loads(manifest.read_text())["comparison_id"] == cid
        assert (STUDY_DIR / f"IOCAS2022_OA_ENERGY_{cid}_rnaseq_report.html").exists()
        assert not (STUDY_DIR / f"IOCAS2022_OA_ENERGY_{cid}_dge_standardized.tsv").exists()


def test_same_lab_gonad_paper_is_not_attached_as_the_source(study):
    assert study["accessions"]["doi"] is None
    assert "cirep" not in study["citation"]


def test_fastq_manifest_covers_every_run_as_paired_end(samplesheet):
    with open(STUDY_DIR / "fastq_manifest.tsv") as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    runs_in_sheet = {r["run_accession"] for r in samplesheet}
    assert {r["run_accession"] for r in manifest} == runs_in_sheet
    assert len(manifest) == 2 * len(runs_in_sheet)
    assert all(len(r["md5"]) == 32 for r in manifest)
