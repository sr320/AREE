"""Curation tests for GEORGE2023_HEATWAVE_PLOIDY.

George et al. 2023 marine heatwave study (PRJNA913164): 72 single-end
QuantSeq 3' Tag-seq libraries, diploid and triploid adult ctenidium, three
treatments. The BioSamples carry no design, so it is joined from the authors'
NCBI submission sheet. These tests pin that join, the prespecified
comparisons and exclusions, the stress_response classification (mortality
was measured per group, not per sequenced oyster), and the Tag-seq settings
the run needs.
"""
from __future__ import annotations

import csv
import json
import re

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "GEORGE2023_HEATWAVE_PLOIDY.yaml"
STUDY_DIR = DATA_DIR / "studies" / "GEORGE2023_HEATWAVE_PLOIDY"
CONFIG = REPO_ROOT / "config" / "GEORGE2023_HEATWAVE_PLOIDY.config"
AUTHORS_EXCLUDED = {"D54", "N56", "X44", "R53", "M43", "N54", "X42", "T62"}
CONTROLS = {"diploid": "diploid_control", "triploid": "triploid_control"}


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


def _tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def samplesheet():
    return _tsv(STUDY_DIR / "samplesheet.tsv")


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_design_is_twelve_per_ploidy_by_treatment_cell(samplesheet):
    prov = json.loads((STUDY_DIR / "ena_provenance.json").read_text())
    assert prov["bioproject"] == "PRJNA913164"
    assert prov["library_layouts"] == ["SINGLE"]
    assert prov["n_runs"] == len(samplesheet) == 72
    assert set(prov["group_sizes"].values()) == {12} and len(prov["group_sizes"]) == 6


def test_design_table_is_derived_from_the_pinned_authors_sheet(samplesheet):
    """design.tsv only splits the authors' description column; check it did
    so faithfully, and that the provenance records what was joined."""
    authors = {r["isolate"]: r["description"] for r in _tsv(STUDY_DIR / "authors_ncbi_attributes.tsv")}
    design = {r["isolate"]: r for r in _tsv(STUDY_DIR / "design.tsv")}
    assert set(authors) == set(design)
    words = {"control": "control", "single stressor": "single_stressor", "multiple stressor": "multiple_stressor"}
    for isolate, desc in authors.items():
        ploidy, treatment = [x.strip() for x in desc.split(",")]
        assert design[isolate]["ploidy"] == ploidy
        assert design[isolate]["treatment"] == words[treatment]

    from common import sha256_file

    recorded = json.loads((STUDY_DIR / "ena_provenance.json").read_text())["design_table"]
    assert recorded["join_attribute"] == "isolate"
    assert recorded["sha256"] == sha256_file(STUDY_DIR / "design.tsv")
    assert "24e7abd1c822504cdff9bbd5dab928f4106b38f9" in recorded["source"]
    for row in samplesheet:
        assert row["condition"] == f"{row['ploidy']}_{row['treatment']}"
        assert row["isolate"] == row["sample_alias"]


def test_each_stressor_is_contrasted_with_its_ploidy_matched_control(study):
    ids = [c["comparison_id"] for c in study["comparisons"]]
    assert ids == ["diploid_30c_vs_control_20c", "triploid_30c_vs_control_20c",
                   "diploid_30c_44c_emersion_vs_control_20c", "triploid_30c_44c_emersion_vs_control_20c"]
    for comp in study["comparisons"]:
        ploidy = comp["comparison_id"].split("_")[0]
        assert comp["control_condition"].startswith(f"{ploidy}, control")
        assert comp["treatment_condition"].startswith(ploidy)
        expected = "multi_stressor" if "emersion" in comp["comparison_id"] else "temperature"
        assert comp["stressor_standardized"] == expected


def test_only_the_prespecified_diploid_30c_comparison_pools(study):
    primary = [c["comparison_id"] for c in study["comparisons"] if c.get("meta_analysis_primary")]
    assert primary == ["diploid_30c_vs_control_20c"]


def test_group_level_mortality_is_not_claimed_as_resilience(study):
    for comp in study["comparisons"]:
        assert comp["phenotype"] == "thermal_tolerance"
        assert comp["resilience_classification"] == "stress_response"
        assert comp["phenotype_direction"] == "not_applicable"
        assert "not per sequenced oyster" in comp["phenotype_unit"]


def test_exclusions_follow_the_authors_and_sample_sizes_add_up(study, samplesheet):
    condition = {r["sample_id"]: r["condition"] for r in samplesheet}
    config_excluded = set(re.search(r"exclude_samples\s*=\s*'([^']*)'", CONFIG.read_text()).group(1).split(","))
    assert config_excluded == AUTHORS_EXCLUDED
    listed = set()
    for comp in study["comparisons"]:
        ploidy = comp["comparison_id"].split("_")[0]
        treatment = "multiple_stressor" if "emersion" in comp["comparison_id"] else "single_stressor"
        arms = {CONTROLS[ploidy], f"{ploidy}_{treatment}"}
        in_arms = {s for s in AUTHORS_EXCLUDED if condition[s] in arms}
        excluded = {e["sample_id"] for e in comp["excluded_samples"]}
        assert excluded == in_arms, comp["comparison_id"]
        assert comp["sample_size"] == 24 - len(excluded)
        listed |= excluded
    assert listed == AUTHORS_EXCLUDED


def test_config_runs_single_end_tag_seq():
    text = CONFIG.read_text()
    assert re.search(r"single_end\s*=\s*true", text)
    assert re.search(r"tag_seq\s*=\s*true", text)


def test_fastq_manifest_covers_every_run_single_end(samplesheet):
    manifest = _tsv(STUDY_DIR / "fastq_manifest.tsv")
    assert {r["run_accession"] for r in manifest} == {r["run_accession"] for r in samplesheet}
    assert len(manifest) == 72
    assert all(len(r["md5"]) == 32 for r in manifest)


def test_replicates_pass_dispersion_qc():
    rows = _tsv(STUDY_DIR / "replicate_dispersion_qc.tsv")
    assert len(rows) == 6
    assert all(r["passes"] == "True" for r in rows)
    assert all(0.02 < float(r["median_dispersion"]) < 0.1 for r in rows)
    # The QC was run on the retained libraries only.
    assert not AUTHORS_EXCLUDED & {s for r in rows for s in r["samples"].split(",")}


def test_each_comparison_was_run_as_registered(study):
    """The manifest's recorded arms must match the comparison it is filed
    under: a config once relabelled every override run as the primary."""
    for comp in study["comparisons"]:
        cid = comp["comparison_id"]
        params = json.loads((STUDY_DIR / f"GEORGE2023_HEATWAVE_PLOIDY_{cid}_workflow_manifest.json").read_text())["parameters"]
        ploidy = cid.split("_")[0]
        treatment = "multiple_stressor" if "emersion" in cid else "single_stressor"
        assert params["comparison_id"] == cid
        assert params["control_levels"] == [CONTROLS[ploidy]]
        assert params["treatment_levels"] == [f"{ploidy}_{treatment}"]
        assert params["single_end"] is True and params["tag_seq"] is True
        assert set(params["exclude_samples"]) == AUTHORS_EXCLUDED
        assert (REPO_ROOT / comp["results_file"]).exists()
        assert (STUDY_DIR / f"GEORGE2023_HEATWAVE_PLOIDY_{cid}_rnaseq_report.html").exists()


def test_standardized_tables_carry_standard_errors(study):
    """Raw reanalysis exists to give poolable evidence: lfcSE and p-values."""
    for comp in study["comparisons"]:
        rows = _tsv(REPO_ROOT / comp["results_file"])
        assert len(rows) > 20000
        assert all(r["lfcSE"] not in ("", "NA") for r in rows[:1000])
        assert sum(r["pvalue"] not in ("", "NA") for r in rows) > 0.95 * len(rows)
