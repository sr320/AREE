"""Curation tests for MA2022_AMMONIA_ZEBRA.

The fourth real Manila clam study and the first with a selected strain: Zebra
II (B) vs ordinary clams (P), unexposed (D) or after a 96 h total-ammonia-
nitrogen challenge (C96), three whole juveniles per cell. These tests pin the
design decoded from the sample aliases, the four prespecified comparisons
(one meta_analysis_primary per contrast type), the direction convention of
the strain contrasts, and the flags that keep the unverified tolerance claim
visible. Tests on the committed run belong here once it exists.
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_ID = "MA2022_AMMONIA_ZEBRA"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
CONFIG_PATH = REPO_ROOT / "config" / f"{STUDY_ID}.config"

# comparison_id -> (control_level, treatment_level, contrast_type, primary)
COMPARISONS = {
    "ordinary_tan96h_vs_ordinary_control": ("P_D", "P_C96", "exposure_response", True),
    "zebra_tan96h_vs_zebra_control": ("B_D", "B_C96", "exposure_response", False),
    "zebra_vs_ordinary_control": ("P_D", "B_D", "baseline_phenotype_difference", True),
    "zebra_vs_ordinary_tan96h": ("P_C96", "B_C96", "phenotype_difference_under_exposure", True),
}


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def samplesheet():
    with open(STUDY_DIR / "samplesheet.tsv", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def ena_provenance():
    return json.loads((STUDY_DIR / "ena_provenance.json").read_text())


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_design_is_decoded_from_the_sample_aliases(samplesheet, ena_provenance):
    """<strain><treatment>_<rep>: B/P = Zebra II/ordinary, D/C96 = control/96 h TAN."""
    assert ena_provenance["bioproject"] == "PRJNA733332"
    assert ena_provenance["run_labels"]["field"] == "sample_alias"
    assert len(samplesheet) == ena_provenance["n_runs"] == ena_provenance["n_samples"] == 12
    assert {r["library_layout"] for r in samplesheet} == {"PAIRED"}
    for row in samplesheet:
        m = re.fullmatch(r"([BP])(D|C96)_([123])", row["sample_alias"])
        assert m, row["sample_alias"]
        assert (row["strain"], row["tan"]) == m.group(1, 2)
        assert row["condition"] == f"{row['strain']}_{row['tan']}"
        assert row["tissue"] == "Complete individual"


def test_three_distinct_clams_per_strain_and_treatment(samplesheet):
    cells = Counter(r["condition"] for r in samplesheet)
    assert cells == dict.fromkeys(("B_D", "B_C96", "P_D", "P_C96"), 3)
    assert len({r["sample_accession"] for r in samplesheet}) == 12
    assert len({r["isolate"] for r in samplesheet}) == 12


def test_four_prespecified_comparisons_one_primary_per_contrast_type(study, samplesheet):
    comps = {c["comparison_id"]: c for c in study["comparisons"]}
    assert set(comps) == set(COMPARISONS)
    levels = {r["condition"] for r in samplesheet}
    primaries = Counter()
    for cid, (ctrl, trt, ctype, primary) in COMPARISONS.items():
        comp = comps[cid]
        assert comp["contrast_type"] == ctype
        assert bool(comp.get("meta_analysis_primary", False)) is primary
        assert f"condition {ctrl})" in comp["control_condition"]
        assert f"condition {trt})" in comp["treatment_condition"]
        assert {ctrl, trt} <= levels
        assert comp["phenotype"] == "toxicant_tolerance"
        assert comp["stressor_standardized"] == "pollutant_exposure"
        assert comp["tissue"] == "whole_animal" and comp["life_stage"] == "juvenile"
        assert comp["sample_size"] == 6 and comp["biological_replicates"] == 3
        primaries[ctype] += primary
    assert all(n == 1 for n in primaries.values())


def test_strain_contrasts_are_resilience_with_zebra_as_treatment(study):
    """Positive log2FC must mean higher in Zebra II; exposure responses claim no tolerance."""
    for comp in study["comparisons"]:
        if comp["contrast_type"] == "exposure_response":
            assert comp["resilience_classification"] == "exposure_only"
            assert comp["phenotype_direction"] == "not_applicable"
        else:
            assert comp["resilience_classification"] == "resilience"
            assert comp["phenotype_direction"] == "increase"
            assert comp["treatment_condition"].startswith("Zebra II")
            assert "not verified" in comp["phenotype_unit"]


def test_unverified_tolerance_and_dose_confound_stay_visible(study):
    assert "ambiguous_phenotype_definition" in study["quality_flags"]
    assert "low_replication" not in study["quality_flags"]  # n = 3 is not < 3
    lim = study["limitations"]
    assert "paywalled" in lim and "different TAN concentrations" in lim
    under = next(c for c in study["comparisons"]
                 if c["contrast_type"] == "phenotype_difference_under_exposure")
    assert "not the same dose" in under["exposure_intensity"]


def test_config_runs_the_full_length_paired_path():
    text = CONFIG_PATH.read_text(encoding="utf-8")
    assert not re.search(r"^\s*tag_seq\s*=\s*true\b", text, re.M)
    assert re.search(r"^\s*single_end\s*=\s*false\b", text, re.M)
    assert re.search(r"^\s*control_level\s*=\s*'P_D'", text, re.M)
    assert re.search(r"^\s*treatment_level\s*=\s*'P_C96'", text, re.M)
    assert "GCF_026571515.1" in text
    block = re.search(r"rnaseq \{(.*?)\}", text, re.S).group(1)
    assert re.search(r"comparison_id\s*=\s*null", block)
