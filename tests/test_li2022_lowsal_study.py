"""Curation tests for LI2022_LOWSAL_VIBRIO.

Li et al. 2022 (PRJNA756710): Vibrio-injected adult gill at 10 vs 30 ppt,
12 h and 48 h, one pooled library per tank. AREE's first real salinity
evidence. Salinity and timepoint are only in the BioSample titles, so these
tests re-derive the design from the verbatim titles kept in design.tsv, and
pin the prespecified comparisons, the tank-level replication, and the
disease (not resilience) classification.
"""
from __future__ import annotations

import csv
import json
import re

import pytest
import yaml

from common import DATA_DIR, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "LI2022_LOWSAL_VIBRIO.yaml"
STUDY_DIR = DATA_DIR / "studies" / "LI2022_LOWSAL_VIBRIO"
TITLE = re.compile(r"post-infection (\d+)h\(replicate=biological replicate-(\d)\)\s+of .* "
                   r"Infected by Vibrio alginolyticus at (low|normal) salini?ty$")


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


def test_design_table_reproduces_the_deposited_titles():
    """Every row's salinity and timepoint follow from its verbatim title
    (the deposit spells 'salinty' in some titles; both spellings parse)."""
    rows = _tsv(STUDY_DIR / "design.tsv")
    assert len(rows) == 14
    for row in rows:
        hours, replicate, salinity = TITLE.search(row["biosample_title"]).groups()
        assert row["timepoint"] == f"{hours}h"
        assert row["title_replicate"] == replicate
        assert row["salinity"] == ("10ppt" if salinity == "low" else "30ppt")


def test_sample_names_agree_with_the_parsed_design(samplesheet):
    """An independent check: C10-n are low salinity, C23-n normal; -1..3 are
    12 h and -4..6 are 48 h; C23hN are the 0 h baselines."""
    for row in samplesheet:
        m = re.fullmatch(r"C(10|23)(-|h)(\d)", row["sample_alias"])
        assert m, row["sample_alias"]
        group, sep, n = m.groups()
        assert row["salinity"] == ("10ppt" if group == "10" else "30ppt")
        expected = "0h" if sep == "h" else ("12h" if int(n) <= 3 else "48h")
        assert row["timepoint"] == expected


def test_design_provenance_records_the_join(samplesheet):
    prov = json.loads((STUDY_DIR / "ena_provenance.json").read_text())
    assert prov["bioproject"] == "PRJNA756710"
    assert prov["group_sizes"] == {"10ppt_12h": 3, "10ppt_48h": 3, "30ppt_0h": 2, "30ppt_12h": 3, "30ppt_48h": 3}
    assert prov["design_table"]["join_attribute"] == "ID"
    assert "TITLE" in prov["design_table"]["source"]


def test_comparisons_are_low_vs_normal_salinity_at_each_timepoint(study):
    ids = [c["comparison_id"] for c in study["comparisons"]]
    assert ids == ["lowsal_10ppt_vs_30ppt_vibrio_48h", "lowsal_10ppt_vs_30ppt_vibrio_12h"]
    for comp in study["comparisons"]:
        assert comp["stressor_standardized"] == "salinity"
        assert comp["sample_size"] == 6 and comp["biological_replicates"] == 3
        assert "3 tanks" in comp["control_condition"] and "3 tanks" in comp["treatment_condition"]


def test_only_48h_pools_because_timepoints_share_tanks(study):
    primary = [c["comparison_id"] for c in study["comparisons"] if c.get("meta_analysis_primary")]
    assert primary == ["lowsal_10ppt_vs_30ppt_vibrio_48h"]


def test_group_mortality_is_disease_evidence_not_resilience(study):
    for comp in study["comparisons"]:
        assert comp["phenotype"] == "disease_susceptibility"
        assert comp["resilience_classification"] == "disease"
        assert comp["phenotype_direction"] == "increase"
        assert "not per sequenced oyster" in comp["phenotype_unit"]
    assert "ambiguous_phenotype_definition" in study["quality_flags"]


def test_fastq_manifest_covers_every_run_paired(samplesheet):
    manifest = _tsv(STUDY_DIR / "fastq_manifest.tsv")
    assert {r["run_accession"] for r in manifest} == {r["run_accession"] for r in samplesheet}
    assert len(manifest) == 28
    assert all(len(r["md5"]) == 32 for r in manifest)
