"""Curation tests for LSU2024_HYPOX_POP.

The first real eastern oyster (Crassostrea virginica) study: gill 3' Tag-seq
of a more hypoxia-tolerant Texas (TX) and a less tolerant Louisiana (LA)
population under hypoxia or control, sampled on days 1-3, 4 oysters per
population x treatment x day. Unpublished deposit (PRJNA1153018). These tests
pin the design decoded from the BioSample isolate attribute, the four
prespecified comparisons and their classes, the day covariate and Tag-seq
path in every committed run, the ENA MD5 check on every FASTQ, the replicate
dispersion gate, and that the stated-not-measured tolerance stays flagged.
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_ID = "LSU2024_HYPOX_POP"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
CONFIG_PATH = REPO_ROOT / "config" / f"{STUDY_ID}.config"
GROUPS = {"LA_control", "LA_hypoxia", "TX_control", "TX_hypoxia"}
LEVELS = {
    "la_hypoxia_vs_control": ("LA_control", "LA_hypoxia"),
    "tx_hypoxia_vs_control": ("TX_control", "TX_hypoxia"),
    "tx_vs_la_hypoxia": ("LA_hypoxia", "TX_hypoxia"),
    "tx_vs_la_control": ("LA_control", "TX_control"),
}
CLASSES = {
    "la_hypoxia_vs_control": ("exposure_response", "stress_response", True),
    "tx_hypoxia_vs_control": ("exposure_response", "stress_response", False),
    "tx_vs_la_hypoxia": ("phenotype_difference_under_exposure", "resilience", True),
    "tx_vs_la_control": ("baseline_phenotype_difference", "resilience", True),
}


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def samplesheet():
    with open(STUDY_DIR / "samplesheet.tsv", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _manifest(comparison_id):
    path = STUDY_DIR / f"{STUDY_ID}_{comparison_id}_workflow_manifest.json"
    return json.loads(path.read_text())


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_design_is_decoded_from_the_isolate_attribute(samplesheet):
    assert len(samplesheet) == 48
    assert len({r["sample_accession"] for r in samplesheet}) == 48
    assert {r["library_layout"] for r in samplesheet} == {"SINGLE"}
    for row in samplesheet:
        m = re.fullmatch(r"(TX|LA) (day[123]) (control|hypoxia)_([1-4])", row["isolate"])
        assert m, row["isolate"]
        assert (row["population"], row["day"], row["treatment"], row["replicate"]) == m.groups()
        assert row["condition"] == f"{m.group(1)}_{m.group(3)}"
        assert row["tissue"] == "gill"


def test_four_oysters_per_population_treatment_and_day(samplesheet):
    cells = Counter((r["condition"], r["day"]) for r in samplesheet)
    assert set(cells.values()) == {4}
    assert {c for c, _ in cells} == GROUPS and len(cells) == 12


def test_comparisons_and_classes_are_as_prespecified(study):
    comps = {c["comparison_id"]: c for c in study["comparisons"]}
    assert set(comps) == set(CLASSES)
    for cid, (contrast, klass, primary) in CLASSES.items():
        comp = comps[cid]
        assert comp["contrast_type"] == contrast
        assert comp["resilience_classification"] == klass
        assert bool(comp.get("meta_analysis_primary", False)) is primary
        assert comp["phenotype"] == "hypoxia_tolerance"
        assert comp["stressor_standardized"] == "hypoxia"
        assert comp["tissue"] == "gill"
        assert comp["sample_size"] == 24 and comp["biological_replicates"] == 12
        assert (REPO_ROOT / comp["results_file"]).exists()


def test_stated_tolerance_stays_flagged(study):
    assert study["species"] == "Crassostrea virginica"
    assert {"ambiguous_phenotype_definition", "family_level_replication"} <= set(
        study["quality_flags"]
    )
    for comp in study["comparisons"]:
        assert comp["exposure_intensity"] is None
    assert study["accessions"]["doi"] is None
    assert study["analysis_status"] == "harmonized"
    assert study["qc_status"] == "passed_with_warnings"


@pytest.mark.parametrize("comparison_id", sorted(LEVELS))
def test_the_committed_runs_are_the_ones_registered(comparison_id):
    """Tag-seq path, day covariate, 12 vs 12 oysters, reference as fetched."""
    manifest = _manifest(comparison_id)
    control, treatment = LEVELS[comparison_id]
    params = manifest["parameters"]
    assert params["tag_seq"] is True and params["single_end"] is True
    assert params["control_levels"] == [control]
    assert params["treatment_levels"] == [treatment]
    design = manifest["deseq2_design"]
    assert design["design"] == "~day + condition"
    assert design["covariates_used"] == ["day"]
    assert design["length_model"].startswith("none")
    assert design["n_libraries_or_units"] == 24
    assert design["residual_df"] == 20
    assert all(f.get("matches_bundle_provenance", True) for f in manifest["reference"]["files"])
    assert manifest["warnings"] == []


def test_replicate_dispersion_qc_passes():
    with open(STUDY_DIR / "replicate_dispersion_qc.tsv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["condition"] for r in rows} == GROUPS
    assert all(r["passes"] == "True" for r in rows)


def test_config_runs_tag_seq_with_day_and_la_exposure_as_default():
    text = CONFIG_PATH.read_text(encoding="utf-8")
    assert re.search(r"^\s*tag_seq\s*=\s*true\b", text, re.M)
    assert re.search(r"^\s*single_end\s*=\s*true\b", text, re.M)
    assert re.search(r"^\s*covariates\s*=\s*'day'", text, re.M)
    assert re.search(r"^\s*control_level\s*=\s*'LA_control'", text, re.M)
    assert re.search(r"^\s*treatment_level\s*=\s*'LA_hypoxia'", text, re.M)
    block = re.search(r"rnaseq \{(.*?)\}", text, re.S).group(1)
    assert re.search(r"comparison_id\s*=\s*null", block)
