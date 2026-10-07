"""Curation tests for LUTIER2022_OA_TIPPING.

A 15-tank pH gradient (one tank per pH, 5 juvenile oysters per tank), not yet
reanalyzed. These tests pin the design read from ENA and the decisions that
are easy to erode later: every run is kept even though 31 are mislabelled
AMPLICON, tank is a condition level so the two pHT 6.5 tanks stay apart, the
prespecified arms, exposure_only classification, and that nothing is pooled
until tank is modelled as the unit of replication.
"""
from __future__ import annotations

import csv
import json
import re

import pytest
import yaml

from common import DATA_DIR, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "LUTIER2022_OA_TIPPING.yaml"
STUDY_DIR = DATA_DIR / "studies" / "LUTIER2022_OA_TIPPING"
CONTROL_LEVELS = {"pHT_7_8_3", "pHT_7_7_5", "pHT_7_6_8"}
TREATMENT_LEVELS = {"pHT_6_7_7", "pHT_6_6_16", "pHT_6_5_9", "pHT_6_5_bis_15"}


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


def _ph(level: str) -> float:
    whole, tenth = re.match(r"pHT_(\d)_(\d)_", level).groups()
    return int(whole) + int(tenth) / 10


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_samplesheet_keeps_the_runs_mislabelled_amplicon(samplesheet, ena_provenance):
    """All 76 runs are RNA-seq; 31 are deposited as AMPLICON. Keep every one."""
    assert ena_provenance["bioproject"] == "PRJNA735889"
    assert ena_provenance["n_runs"] == ena_provenance["n_samples"] == len(samplesheet) == 76
    assert set(ena_provenance["library_strategies"]) == {"AMPLICON", "RNA-Seq"}
    assert len({r["sample_alias"] for r in samplesheet}) == 76


def test_every_condition_level_is_one_tank(samplesheet):
    tanks_per_level: dict[str, set[str]] = {}
    for row in samplesheet:
        tanks_per_level.setdefault(row["condition"], set()).add(row["tank"])
    assert len(tanks_per_level) == 15
    assert all(len(t) == 1 for t in tanks_per_level.values())
    # The two pHT 6.5 tanks (nominal 6.4 and 6.5) must not be merged.
    assert {"pHT_6_5_9", "pHT_6_5_bis_15"} <= set(tanks_per_level)


def test_prespecified_arms_bracket_the_tipping_points(study, samplesheet):
    comps = study["comparisons"]
    assert [c["comparison_id"] for c in comps] == ["ph_le6_7_vs_ph_ge7_6_23d"]
    comp = comps[0]
    assert comp["meta_analysis_primary"] is True

    levels = {r["condition"] for r in samplesheet}
    assert levels >= CONTROL_LEVELS | TREATMENT_LEVELS
    for level in CONTROL_LEVELS | TREATMENT_LEVELS:
        assert level in comp["control_condition"] + comp["treatment_condition"]
    # Paper tipping points: 7.3-6.9. Arms sit outside that band on both sides.
    assert min(_ph(lv) for lv in CONTROL_LEVELS) >= 7.6
    assert max(_ph(lv) for lv in TREATMENT_LEVELS) <= 6.7

    in_arms = [r for r in samplesheet if r["condition"] in CONTROL_LEVELS | TREATMENT_LEVELS]
    assert len(in_arms) == comp["sample_size"] == 35
    # Replication is counted in tanks, not oysters.
    assert comp["biological_replicates"] == min(len(CONTROL_LEVELS), len(TREATMENT_LEVELS)) == 3


def test_exposure_contrast_is_exposure_only(study):
    comp = study["comparisons"][0]
    assert comp["stressor_standardized"] == "ocean_acidification"
    assert comp["phenotype"] == "acidification_tolerance"
    assert comp["resilience_classification"] == "exposure_only"
    assert comp["phenotype_direction"] == "not_applicable"
    assert comp["tissue"] == "whole_animal"
    assert comp["life_stage"] == "juvenile"


def test_publication_is_the_confirmed_one(study):
    assert study["accessions"]["doi"] == "10.1111/gcb.16101"
    # The 2025 J Exp Biol paper reuses these data; it is not this study's DOI.
    assert "10.1242/jeb.249458" not in study["citation"]


def test_nothing_runs_or_pools_until_tank_is_the_unit(study):
    assert study["qc_status"] == "not_started"
    assert study["analysis_status"] == "not_started"
    assert all(c["results_file"] is None for c in study["comparisons"])
    assert "tank" in study["limitations"].lower()
