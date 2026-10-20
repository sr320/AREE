"""Curation tests for LI2024_VIBRIO_GILL.

C. gigas gill at 0, 6 and 48 h after Vibrio alginolyticus infection, n=3, one
BioSample per oyster. These tests pin the design read from ENA and the
decisions that are easy to erode later: 0h_replicate3 is excluded because it is
mantle, not gill, on evidence independent of the infection contrast; the faint
mantle signal in 48h_replicate3 does not justify excluding it; the only control
is the 0 h baseline, so both comparisons are disease exposure contrasts; and
the shell-matrix genes no longer read as "down at 48 h" once the mantle
baseline is gone.
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_ID = "LI2024_VIBRIO_GILL"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
CONFIG_PATH = REPO_ROOT / "config" / f"{STUDY_ID}.config"
TISSUE_DIR = STUDY_DIR / "tissue_identity_check"
EXCLUDED = "0h_replicate3"
COMPARISONS = {"valginolyticus_48h_vs_0h": "48h", "valginolyticus_6h_vs_0h": "6h"}
SHELL_GENES = ("105319088", "105346498", "105338611", "105331575")


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def comparisons(study):
    return {c["comparison_id"]: c for c in study["comparisons"]}


@pytest.fixture(scope="module")
def samplesheet():
    with open(STUDY_DIR / "samplesheet.tsv", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def _results(comp):
    return {r["gene_id"]: r for r in _read_tsv(REPO_ROOT / comp["results_file"])}


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_three_distinct_gill_oysters_per_timepoint(samplesheet):
    assert len(samplesheet) == 9
    assert len({r["sample_accession"] for r in samplesheet}) == 9
    assert Counter(r["condition"] for r in samplesheet) == {"0h": 3, "6h": 3, "48h": 3}
    assert {r["tissue"] for r in samplesheet} == {"gill"}  # as deposited; see exclusion
    assert {r["library_layout"] for r in samplesheet} == {"PAIRED"}


def test_two_comparisons_against_the_0h_baseline(comparisons):
    """No time-matched uninfected oysters exist: infection is confounded with time,
    so neither comparison can be a resistance phenotype."""
    assert set(comparisons) == set(COMPARISONS)
    primary = [cid for cid, c in comparisons.items() if c.get("meta_analysis_primary")]
    assert primary == ["valginolyticus_48h_vs_0h"]
    for c in comparisons.values():
        assert c["resilience_classification"] == "disease"
        assert c["contrast_type"] == "exposure_response"
        assert c["phenotype"] == "disease_susceptibility"
        assert c["stressor_standardized"] == "pathogen_challenge"
        assert c["tissue"] == "gill"
        assert c["sample_size"] == 5
        assert c["biological_replicates"] == 2


def test_mantle_baseline_is_excluded_from_both_comparisons(comparisons):
    for c in comparisons.values():
        excluded = c["excluded_samples"]
        assert [e["sample_id"] for e in excluded] == [EXCLUDED]
        assert "tissue_identity_check" in excluded[0]["evidence"]
        assert "not the infection contrast" in excluded[0]["evidence"]
    text = CONFIG_PATH.read_text(encoding="utf-8")
    assert re.search(rf"^\s*exclude_samples\s*=\s*'{EXCLUDED}'", text, re.M)


def test_exclusion_rests_on_tissue_markers_not_on_the_result():
    """The excluded library expresses shell-matrix genes at mantle levels and lacks
    gill prestin. 48h_replicate3 carries about 1% of that signal and is kept."""
    rows = {r["gene_id"]: r for r in _read_tsv(TISSUE_DIR / "marker_genes_tpm.tsv")}
    others = [s for s in rows[SHELL_GENES[0]] if s.endswith(tuple("123")) and s != EXCLUDED]
    assert len(others) == 8
    for gene in SHELL_GENES:
        assert float(rows[gene][EXCLUDED]) > 1_000
        assert max(float(rows[gene][s]) for s in others) < 0.05 * float(rows[gene][EXCLUDED])
    prestin = rows["105336253"]
    assert float(prestin[EXCLUDED]) < min(float(prestin[s]) for s in others)

    pca = {r["sample_id"]: r for r in _read_tsv(TISSUE_DIR / "pca_and_correlation.tsv")}
    assert float(pca[EXCLUDED]["PC1"]) > 100
    assert all(float(pca[EXCLUDED][k]) < 0 for k in ("mean_r_0h", "mean_r_6h", "mean_r_48h"))
    assert float(pca["48h_replicate3"]["mean_r_48h"]) > 0.75


@pytest.mark.parametrize("cid,level", COMPARISONS.items())
def test_manifest_records_the_exclusion(cid, level):
    manifest = json.loads((STUDY_DIR / f"{STUDY_ID}_{cid}_workflow_manifest.json").read_text())
    params = manifest["parameters"]
    assert params["exclude_samples"] == [EXCLUDED]
    assert params["control_levels"] == ["0h"]
    assert params["treatment_levels"] == [level]
    assert params["single_end"] is False
    assert manifest["deseq2_design"]["n_libraries_or_units"] == 5
    assert manifest["deseq2_design"]["residual_df"] == 3
    assert all(f.get("matches_bundle_provenance", True) for f in manifest["reference"]["files"])
    assert manifest["warnings"] == [
        f"samples_excluded_by_curation: {EXCLUDED} (see the study registry for the reason)"]


@pytest.mark.parametrize("cid,expected", [
    ("valginolyticus_48h_vs_0h", (26_652, 3_753, 1_459)),
    ("valginolyticus_6h_vs_0h", (25_888, 1_133, 430)),
])
def test_results_carry_standard_errors(comparisons, cid, expected):
    results = _results(comparisons[cid])
    n, n_sig, n_up = expected
    assert len(results) == n
    assert all(r["lfcSE"] not in ("", "NA") for r in results.values())
    sig = [r for r in results.values() if r["padj"] not in ("", "NA") and float(r["padj"]) < 0.05]
    assert len(sig) == n_sig
    assert sum(float(r["log2FoldChange"]) > 0 for r in sig) == n_up


def test_shell_genes_are_not_an_infection_signal_at_48h(comparisons):
    """With the mantle baseline included they read about -7.6 at 48 h."""
    results = _results(comparisons["valginolyticus_48h_vs_0h"])
    for gene in SHELL_GENES:
        assert abs(float(results[gene]["log2FoldChange"])) < 1
        assert float(results[gene]["padj"]) > 0.5


def test_replicate_dispersion_passes_after_the_exclusion():
    rows = {r["condition"]: r for r in _read_tsv(STUDY_DIR / "replicate_dispersion_qc.tsv")}
    assert set(rows) == {"0h", "6h", "48h"}
    assert EXCLUDED not in rows["0h"]["samples"]
    assert int(rows["0h"]["n_samples"]) == 2
    assert all(r["passes"] == "True" for r in rows.values())


def test_caveats_stay_on_record(study):
    assert "low_replication" in study["quality_flags"]
    assert "confounded" in study["limitations"]
    assert "inferred" in study["citation"]
    assert study["qc_status"] == "passed_with_warnings"
    assert study["analysis_status"] == "harmonized"
