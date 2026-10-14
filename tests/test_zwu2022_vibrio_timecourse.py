"""Curation tests for ZWU2022_VIBRIO_TIMECOURSE.

Vibrio exposure time course (PRJNA852380), juvenile digestive gland, untreated
controls and 1-24 h, three libraries each; no publication. These tests pin
the design read from the deposit, the pre-download genotype screen, and the
prespecified comparisons.
"""
from __future__ import annotations

import csv
import json

import pytest
import yaml

from common import DATA_DIR, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "ZWU2022_VIBRIO_TIMECOURSE.yaml"
STUDY_DIR = DATA_DIR / "studies" / "ZWU2022_VIBRIO_TIMECOURSE"
TIMES = (24, 12, 6, 3, 1)


def _tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_design_matches_biosample_isolate():
    rows = _tsv(STUDY_DIR / "samplesheet.tsv")
    prov = json.loads((STUDY_DIR / "ena_provenance.json").read_text())
    assert prov["n_runs"] == len(rows) == 18
    assert set(prov["group_sizes"].values()) == {3}
    for row in rows:
        expected = "normal" if row["condition"] == "control" else f"vibrio stress {row['condition'][1:]}"
        assert row["isolate"] == expected


def test_genotype_screen_shows_distinct_animals():
    text = (STUDY_DIR / "replicate_genotype_check.txt").read_text()
    within = [line for line in text.splitlines() if line.strip().startswith("within")][0]
    observed, noise = float(within.split("observed")[1].split()[0]), float(within.split("noise")[1].split()[0])
    assert observed > 2 * noise


def test_each_time_point_vs_control_with_24h_primary(study):
    ids = [c["comparison_id"] for c in study["comparisons"]]
    assert ids == [f"vibrio_{t}h_vs_control" for t in TIMES]
    primary = [c["comparison_id"] for c in study["comparisons"] if c.get("meta_analysis_primary")]
    assert primary == ["vibrio_24h_vs_control"]
    for comp in study["comparisons"]:
        assert comp["resilience_classification"] == "disease"
        assert comp["exposure_intensity"] is None  # not deposited, no paper
