"""Curation tests for UW2016_HEATSHOCK_MANTLE.

Unpublished Roberts-lab data: mantle from three adult oysters before (15 C) and
after a 1 h heat shock at 40 C. These tests pin the design read from ENA and
the decisions that are easy to erode later: each oyster is its own control
(~ oyster + condition, 2 residual df), the single prespecified comparison is
exposure_only, the low-mapping 6M-HS library is kept, and the inducible HSP70
B2 rises, which is the check that the C/HS labels are not swapped.
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_ID = "UW2016_HEATSHOCK_MANTLE"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
CONFIG_PATH = REPO_ROOT / "config" / f"{STUDY_ID}.config"
COMPARISON_ID = "heat_shock_40c_1h_vs_control"
MANIFEST_PATH = STUDY_DIR / f"{STUDY_ID}_{COMPARISON_ID}_workflow_manifest.json"
HSP70_B2 = "105334234"  # LOC105334234, inducible HSP70 B2


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def samplesheet():
    with open(STUDY_DIR / "samplesheet.tsv", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def manifest():
    return json.loads(MANIFEST_PATH.read_text())


@pytest.fixture(scope="module")
def results(study):
    path = REPO_ROOT / study["comparisons"][0]["results_file"]
    with open(path, encoding="utf-8") as fh:
        return {r["gene_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_three_oysters_each_sampled_before_and_after(samplesheet):
    """Six runs, six BioSamples, but only three animals: the pairing is the design."""
    assert len(samplesheet) == 6
    assert len({r["sample_accession"] for r in samplesheet}) == 6
    assert Counter(r["condition"] for r in samplesheet) == {"C": 3, "HS": 3}
    pairs = Counter((r["oyster"], r["condition"]) for r in samplesheet)
    assert set(pairs) == {(o, c) for o in ("2M", "4M", "6M") for c in ("C", "HS")}
    assert all(n == 1 for n in pairs.values())
    assert {r["tissue"] for r in samplesheet} == {"mantle"}
    assert {r["library_layout"] for r in samplesheet} == {"SINGLE"}


def test_config_fits_the_paired_design(samplesheet):
    text = CONFIG_PATH.read_text(encoding="utf-8")
    assert re.search(r"^\s*covariates\s*=\s*'oyster'", text, re.M)
    assert re.search(r"^\s*single_end\s*=\s*true\b", text, re.M)
    levels = set(re.findall(r"(?:control_level|treatment_level)\s*=\s*'([^']+)'", text))
    assert levels == {"C", "HS"} <= {r["condition"] for r in samplesheet}


def test_command_line_comparison_id_is_not_shadowed():
    """The workflow reads rnaseq.comparison_id before comparison_id; set inside the
    rnaseq block it would silently override the top-level one."""
    text = CONFIG_PATH.read_text(encoding="utf-8")
    block = re.search(r"rnaseq \{(.*?)\}", text, re.S).group(1)
    assert re.search(r"comparison_id\s*=\s*null", block)
    assert re.search(rf"^\s*comparison_id\s*=\s*'{COMPARISON_ID}'", text, re.M)


def test_one_primary_comparison_classified_as_exposure(study):
    """Heat shock vs the same animal a week earlier is an exposure. No oyster was
    scored for tolerance, so this must not become resilience evidence."""
    assert [c["comparison_id"] for c in study["comparisons"]] == [COMPARISON_ID]
    comp = study["comparisons"][0]
    assert comp["meta_analysis_primary"] is True
    assert comp["resilience_classification"] == "exposure_only"
    assert comp["contrast_type"] == "exposure_response"
    assert comp["phenotype"] == "thermal_tolerance"
    assert comp["phenotype_direction"] == "not_applicable"
    assert comp["stressor_standardized"] == "temperature"
    assert comp["tissue"] == "mantle"
    assert comp["biological_replicates"] == 3


def test_caveats_stay_on_record(study):
    assert {"low_replication", "single_timepoint"} <= set(study["quality_flags"])
    limitations = study["limitations"]
    assert "2 residual degrees" in limitations
    assert "confounded with time and" in limitations
    assert "2013-03-19" in limitations
    assert study["accessions"]["doi"] is None
    assert study["citation"].startswith("Unpublished")
    assert study["qc_status"] == "passed_with_warnings"
    assert study["analysis_status"] == "harmonized"


def test_manifest_records_the_paired_run(manifest):
    assert manifest["parameters"]["covariates"] == ["oyster"]
    assert manifest["parameters"]["exclude_samples"] == []  # 6M-HS kept, see YAML
    assert manifest["parameters"]["single_end"] is True
    assert manifest["deseq2_design"]["design"] == "~oyster + condition"
    assert manifest["deseq2_design"]["n_libraries_or_units"] == 6
    assert manifest["deseq2_design"]["residual_df"] == 2
    assert manifest["reference"]["genome_assembly"] == "GCF_963853765.1"
    assert all(f.get("matches_bundle_provenance", True) for f in manifest["reference"]["files"])
    assert manifest["warnings"] == []


def test_results_carry_standard_errors(results):
    """Raw reanalysis is the point: every record must be poolable."""
    assert len(results) == 26_555
    assert all(r["lfcSE"] not in ("", "NA") for r in results.values())
    sig = [r for r in results.values() if r["padj"] not in ("", "NA") and float(r["padj"]) < 0.05]
    assert len(sig) == 5_465
    assert sum(float(r["log2FoldChange"]) > 0 for r in sig) == 2_329


def test_inducible_hsp70_rises(results):
    """The label check: if C and HS were swapped, this would fall."""
    hsp = results[HSP70_B2]
    assert float(hsp["log2FoldChange"]) > 3
    assert float(hsp["padj"]) < 0.05


def test_replicate_dispersion_passes_before_and_after_full_depth():
    for name in ("prescreen_replicates.tsv", "replicate_dispersion_qc.tsv"):
        with open(STUDY_DIR / name, encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert {r["condition"] for r in rows} == {"C", "HS"}
        assert all(r["passes"] == "True" for r in rows)
