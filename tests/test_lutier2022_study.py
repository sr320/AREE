"""Curation tests for LUTIER2022_OA_TIPPING.

A 15-tank pH gradient (one tank per pH, 5 juvenile oysters per tank),
reanalyzed with tank as the unit of replication. These tests pin the design
read from ENA and the decisions that are easy to erode later: every run is
kept even though 31 are mislabelled AMPLICON, tank is a condition level so the
two pHT 6.5 tanks stay apart, the prespecified arms, exposure_only
classification, and that the result was fitted on tanks, not oysters.
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


def test_config_runs_the_prespecified_arms_on_tanks():
    from common import REPO_ROOT

    config = (REPO_ROOT / "config" / "LUTIER2022_OA_TIPPING.config").read_text()
    levels = dict(re.findall(r"^\s*(control_level|treatment_level|replicate_unit)\s*=\s*'([^']*)'", config, re.M))
    assert set(levels["control_level"].split(",")) == CONTROL_LEVELS
    assert set(levels["treatment_level"].split(",")) == TREATMENT_LEVELS
    assert levels["replicate_unit"] == "tank"


def test_result_was_fitted_on_tanks_not_oysters(study):
    """The manifest must say the lfcSE came from 7 tanks, not 35 oysters."""
    comp = study["comparisons"][0]
    assert "tank_level_replication" in study["quality_flags"]
    stem = STUDY_DIR / "LUTIER2022_OA_TIPPING_ph_le6_7_vs_ph_ge7_6_23d"
    params = json.loads(stem.with_name(stem.name + "_workflow_manifest.json").read_text())["parameters"]
    assert params["replicate_unit"] == "tank"
    assert set(params["control_levels"]) == CONTROL_LEVELS
    assert set(params["treatment_levels"]) == TREATMENT_LEVELS
    assert comp["results_file"].endswith(stem.name + "_dge_standardized.tsv")


def test_oysters_within_each_tank_are_biological_replicates():
    """Unlike IOCAS2022_OA_ENERGY, every tank's oysters vary like real animals."""
    with open(STUDY_DIR / "replicate_dispersion_qc.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["condition"] for r in rows} == CONTROL_LEVELS | TREATMENT_LEVELS
    assert all(r["passes"] == "True" for r in rows)
    assert all(0.02 < float(r["median_dispersion"]) < 0.1 for r in rows)


def test_standardized_result_has_tank_level_standard_errors(study):
    import pandas as pd

    from common import REPO_ROOT

    df = pd.read_csv(REPO_ROOT / study["comparisons"][0]["results_file"], sep="\t")
    assert {"gene_id", "log2FoldChange", "lfcSE", "pvalue", "padj"} <= set(df.columns)
    # 1,343 genes at padj < 0.05 on tanks; fitting the 35 oysters gave 4,922.
    assert (df["padj"] < 0.05).sum() == 1343
    assert study["qc_status"] == "passed_with_warnings"
