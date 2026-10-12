"""Curation tests for BERNARDINI2024_TIRE.

Bernardini et al. 2024 (PRJNA856813): tire microparticles and leachate,
adult gill on day 41, one oyster per library. AREE's first real pollutant
study. The design comes from run-level experiment titles because the
BioSamples are shared 16S records; these tests pin that parse, that only the
NovaSeq RNA-seq runs enter the comparisons, the 3' tag settings, and the
stress_response classification.
"""
from __future__ import annotations

import csv
import json
import re

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "BERNARDINI2024_TIRE.yaml"
STUDY_DIR = DATA_DIR / "studies" / "BERNARDINI2024_TIRE"
CONFIG = REPO_ROOT / "config" / "BERNARDINI2024_TIRE.config"
CONTROL = "NovaSeq_6000_D41_GI_CTL"
ARMS = {
    "tire_particles_high_vs_ctl_d41_gill": "NovaSeq_6000_D41_GI_MR_H",
    "tire_particles_low_vs_ctl_d41_gill": "NovaSeq_6000_D41_GI_MR_L",
    "tire_leachate_high_vs_ctl_d41_gill": "NovaSeq_6000_D41_GI_LEA_H",
    "tire_leachate_low_vs_ctl_d41_gill": "NovaSeq_6000_D41_GI_LEA_L",
}


def _tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def samplesheet():
    return _tsv(STUDY_DIR / "samplesheet.tsv")


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_every_run_parses_and_the_gill_groups_are_six_novaseq_oysters(samplesheet):
    prov = json.loads((STUDY_DIR / "ena_provenance.json").read_text())
    assert prov["bioproject"] == "PRJNA856813" and prov["n_runs"] == len(samplesheet) == 213
    for level in [CONTROL, *ARMS.values()]:
        rows = [r for r in samplesheet if r["condition"] == level]
        assert len(rows) == 6, level
        assert all(r["sample_id"].startswith("NovaSeq 6000 sequencing: D41_GI_") for r in rows)


def test_16s_runs_never_share_a_condition_with_rnaseq(samplesheet):
    """Each oyster's 16S and RNA-seq runs share a BioSample; the platform in
    the condition keeps the amplicon runs out of every arm."""
    used = {CONTROL, *ARMS.values()}
    for row in samplesheet:
        if row["condition"] in used:
            assert row["condition"].startswith("NovaSeq_6000_")
    miseq = [r for r in samplesheet if r["condition"].startswith("MiSeq_")]
    assert miseq and not ({r["condition"] for r in miseq} & used)


def test_config_runs_single_end_tag_seq_against_the_controls():
    text = CONFIG.read_text()
    assert re.search(r"single_end\s*=\s*true", text)
    assert re.search(r"tag_seq\s*=\s*true", text)
    assert f"control_level   = '{CONTROL}'" in text


def test_comparisons_classification_and_primary(study):
    ids = [c["comparison_id"] for c in study["comparisons"]]
    assert ids == list(ARMS)
    primary = [c["comparison_id"] for c in study["comparisons"] if c.get("meta_analysis_primary")]
    assert primary == ["tire_particles_high_vs_ctl_d41_gill"]
    for comp in study["comparisons"]:
        assert comp["stressor_standardized"] == "pollutant_exposure"
        assert comp["phenotype"] == "toxicant_tolerance"
        assert comp["resilience_classification"] == "stress_response"
        assert comp["exposure_intensity"] is None  # concentrations unknown: paper not read in full
        assert comp["sample_size"] == 12
