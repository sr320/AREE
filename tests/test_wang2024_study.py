"""Curation tests for WANG2024_HEAT_GILL.

C. gigas adult gill, 12 h heat vs no heat, deposited as n=3 per arm with
paired ATAC-seq and a C. angulata arm. It failed QC on pseudo-replication:
replicate dispersion was far too low, and a genotype check showed that each
arm's three libraries share one RNA source. These tests pin the design read
from ENA and keep the failed study out of every pool.
"""
from __future__ import annotations

import csv
import json

import pytest
import yaml

from common import DATA_DIR, STUDIES_DIR

STUDY_ID = "WANG2024_HEAT_GILL"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
CID = "gigas_heat_12h_vs_control"


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def samplesheet():
    with open(STUDY_DIR / "samplesheet.tsv") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _tsv(name: str) -> list[dict]:
    with open(STUDY_DIR / name) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_design_is_parsed_from_library_name(samplesheet):
    """sample_alias is shared by an oyster's RNA-seq and ATAC-seq runs; only
    library_name identifies a library."""
    provenance = json.loads((STUDY_DIR / "ena_provenance.json").read_text())
    assert provenance["bioproject"] == "PRJNA970954"
    assert provenance["n_runs"] == len(samplesheet) == 24
    assert provenance["n_samples"] == 12
    assert provenance["group_sizes"] == dict.fromkeys(
        [f"{a}_{s}_{h}" for a in ("ATAC_Seq", "RNA_Seq") for s in ("angulata", "gigas") for h in ("HC", "NC")], 3)
    assert len({r["sample_id"] for r in samplesheet}) == 24


def test_registered_arm_is_gigas_rnaseq_only(study, samplesheet):
    (comp,) = study["comparisons"]
    assert comp["comparison_id"] == CID
    arm = [r for r in samplesheet if r["condition"] in ("RNA_Seq_gigas_HC", "RNA_Seq_gigas_NC")]
    assert len(arm) == comp["sample_size"] == 6
    assert len({r["sample_accession"] for r in arm}) == 6
    assert comp["resilience_classification"] == "exposure_only"
    assert comp["exposure_duration"] == "12 h"
    assert comp["exposure_intensity"] is None   # temperature not yet read from the methods


def test_failed_qc_keeps_the_study_out_of_every_pool(study):
    assert study["qc_status"] == "failed"
    assert study["analysis_status"] == "in_progress"
    assert not any(c.get("meta_analysis_primary") for c in study["comparisons"])
    assert all(c["results_file"] is None for c in study["comparisons"])
    # The run is archived as evidence, but its standardized table is not committed.
    assert (STUDY_DIR / f"{STUDY_ID}_{CID}_workflow_manifest.json").exists()
    assert not (STUDY_DIR / f"{STUDY_ID}_{CID}_dge_standardized.tsv").exists()


def test_dispersion_qc_records_the_failure():
    rows = {r["condition"]: r for r in _tsv("replicate_dispersion_qc.tsv")}
    assert set(rows) == {"RNA_Seq_gigas_HC", "RNA_Seq_gigas_NC"}
    assert all(r["passes"] == "False" for r in rows.values())
    assert all(float(r["median_dispersion"]) < 0.005 for r in rows.values())


def test_genotype_check_shows_one_source_per_arm():
    """Within an arm, libraries differ about as little as two libraries of one
    RNA would; between arms they differ several times more."""
    rows = _tsv("genotype_check/pairwise_allele_fractions.tsv")
    assert len(rows) == 15
    assert {int(r["n_sites"]) for r in rows} == {9355}
    within = [float(r["ratio_observed_to_same_RNA"]) for r in rows if r["same_arm"] == "True"]
    between = [float(r["ratio_observed_to_same_RNA"]) for r in rows if r["same_arm"] == "False"]
    assert len(within) == 6 and len(between) == 9
    assert max(within) < 2.5 < 5 < min(between)
    assert (STUDY_DIR / "genotype_check" / "run_genocheck.sh").exists()
