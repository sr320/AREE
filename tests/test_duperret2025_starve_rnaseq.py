"""Curation tests for DUPERRET2025_STARVE_RNASEQ.

RNA-seq arm of Duperret et al. 2025 (iScience), deposited in ENA as two
projects split by diet (PRJEB86618 starved, PRJEB86525 fed) and registered as
one study before any reanalysis. These tests pin the merge of the two
fetch-samplesheet outputs, the design parsed from experiment_title, and the
curation decisions that are easy to erode later: both projects listed as
accessions, tank-level replication, the diet contrast at 0 hpc as primary,
and the infection contrasts kept as `disease` evidence.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_ID = "DUPERRET2025_STARVE_RNASEQ"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
PROJECTS = {"PRJEB86618": "5", "PRJEB86525": "3"}  # project -> diet code (5 starved, 3 fed)
TIMES = ("0", "3", "6", "12", "24", "48", "72")


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def samplesheet():
    with open(STUDY_DIR / "samplesheet.tsv") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def provenance():
    return json.loads((STUDY_DIR / "ena_provenance.json").read_text())


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_both_projects_are_merged_and_each_holds_one_diet(samplesheet, provenance):
    assert provenance["bioprojects"] == list(PROJECTS)
    assert provenance["n_runs"] == len(samplesheet) == 83
    per_project = provenance["per_project_provenance"]
    assert set(per_project) == set(PROJECTS)
    run_to_project = {}
    for proj, prov in per_project.items():
        for run in prov["run_labels"]["labels"]:
            run_to_project[run] = proj
    for r in samplesheet:
        assert r["diet_code"] == PROJECTS[run_to_project[r["run_accession"]]]
        assert r["experimental_factor"].split(",")[0] == {"3": "control", "5": "starvation"}[r["diet_code"]]


def test_design_is_families_by_diets_by_seven_times(samplesheet, provenance):
    cells = Counter((r["family"], r["diet_code"], r["hpc"]) for r in samplesheet)
    expected = {(f, d, t) for f in ("F11N", "F14R") for d in ("3", "5") for t in TIMES}
    assert set(cells) == expected
    # One planned library (F14R starved 6 hpc) was never deposited.
    assert {k: v for k, v in cells.items() if v != 3} == {("F14R", "5", "6"): 2}
    assert provenance["library_layouts"] == ["SINGLE"]
    for r in samplesheet:
        assert r["condition"] == f"{r['diet_code']}_{r['hpc']}"
        assert r["sample_id"] == f"P_{r['family']}_{r['diet_code']}_{r['hpc']}_{r['replicate']}"


def test_fastq_manifest_is_single_end_with_md5s(samplesheet):
    with open(STUDY_DIR / "fastq_manifest.tsv") as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["run_accession"] for r in manifest} == {r["run_accession"] for r in samplesheet}
    assert len(manifest) == len(samplesheet)
    assert all(len(r["md5"]) == 32 for r in manifest)


def test_author_excluded_outlier_is_outside_registered_comparisons(study, provenance, samplesheet):
    out = provenance["author_excluded"]
    assert out["sample_id"] == "P_F11N_3_48_2"
    row = next(r for r in samplesheet if r["sample_id"] == out["sample_id"])
    used = {"3_0", "5_0", "3_24", "5_24"}
    assert row["condition"] not in used


def test_accessions_name_both_projects(study):
    accs = study["accessions"]
    assert accs["bioproject"] == "PRJEB86618"
    assert "PRJEB86618" in accs["ena"] and "PRJEB86525" in accs["ena"]
    assert accs["doi"] == "10.1016/j.isci.2025.114333"
    assert "MTBLS13313" not in " ".join(str(v) for v in accs.values() if v)


def test_comparisons_are_tank_replicated_and_classified(study):
    comps = {c["comparison_id"]: c for c in study["comparisons"]}
    assert set(comps) == {"starved_vs_fed_0hpc", "fed_oshv1_24hpc_vs_0hpc", "starved_oshv1_24hpc_vs_0hpc"}
    assert [cid for cid, c in comps.items() if c.get("meta_analysis_primary")] == ["starved_vs_fed_0hpc"]
    diet = comps["starved_vs_fed_0hpc"]
    assert diet["stressor_standardized"] == "nutritional_limitation"
    assert diet["phenotype"] == "disease_resistance"
    assert diet["resilience_classification"] == "resilience"
    for cid in ("fed_oshv1_24hpc_vs_0hpc", "starved_oshv1_24hpc_vs_0hpc"):
        assert comps[cid]["resilience_classification"] == "disease"
        assert comps[cid]["phenotype"] == "disease_susceptibility"
    for c in comps.values():
        assert c["sample_size"] == 12 and c["biological_replicates"] == 3
    assert "tank_level_replication" in study["quality_flags"]


def test_every_comparison_fits_three_tank_units_per_arm(samplesheet):
    # Both families share each tank, so each comparison is fitted on tank
    # units (families summed): 3 vs 3, never 6 vs 6 libraries.
    contrasts = {
        "starved_vs_fed_0hpc": (("3_0", "5_0"), "tank"),
        "fed_oshv1_24hpc_vs_0hpc": (("3_0", "3_24"), "tank_hpc"),
        "starved_oshv1_24hpc_vs_0hpc": (("5_0", "5_24"), "tank_hpc"),
    }
    for (control, treatment), unit in contrasts.values():
        rows = [r for r in samplesheet if r["condition"] in (control, treatment)]
        units = {}
        for r in rows:
            units.setdefault(r[unit], set()).add(r["condition"])
            assert r["tank"] == f"{r['diet']}_t{r['replicate']}"
            assert r["tank_hpc"] == f"{r['tank']}_{r['hpc']}hpc"
        # A unit never spans both arms (the workflow would refuse it).
        assert all(len(arms) == 1 for arms in units.values())
        per_arm = Counter(next(iter(arms)) for arms in units.values())
        assert per_arm == {control: 3, treatment: 3}
        families = Counter(r[unit] for r in rows)
        assert set(families.values()) == {2}


def test_config_is_single_end_and_fits_primary_on_tanks():
    config = (REPO_ROOT / "config" / f"{STUDY_ID}.config").read_text()
    assert "replicate_unit  = 'tank'" in config
    assert "single_end = true" in config
    assert "covariates" not in config.split("params {", 1)[1]


def test_outputs_match_provenance_checksums(provenance):
    for rel, digest in provenance["outputs"].items():
        assert hashlib.sha256((REPO_ROOT / rel).read_bytes()).hexdigest() == digest


# ---- 2026-10-12 reanalysis outputs -------------------------------------------

RUN = {
    "starved_vs_fed_0hpc": (["3_0"], ["5_0"], "tank", False),
    "fed_oshv1_24hpc_vs_0hpc": (["3_0"], ["3_24"], "tank_hpc", True),
    "starved_oshv1_24hpc_vs_0hpc": (["5_0"], ["5_24"], "tank_hpc", True),
}


def test_registry_points_at_committed_results(study):
    assert study["qc_status"] == "passed_with_warnings"
    assert study["analysis_status"] == "harmonized"
    for comp in study["comparisons"]:
        assert comp["results_file"] == (
            f"data/studies/{STUDY_ID}/{STUDY_ID}_{comp['comparison_id']}_dge_standardized.tsv")
        assert (REPO_ROOT / comp["results_file"]).exists()


@pytest.mark.parametrize("comparison_id", list(RUN))
def test_manifest_records_the_registered_fit(comparison_id):
    control, treatment, unit, viral = RUN[comparison_id]
    manifest = json.loads((STUDY_DIR / f"{STUDY_ID}_{comparison_id}_workflow_manifest.json").read_text())
    params = manifest["parameters"]
    assert manifest["mode"] == "raw_reanalysis"
    assert params["control_levels"] == control and params["treatment_levels"] == treatment
    assert params["replicate_unit"] == unit
    assert params["single_end"] is True
    assert params["se_adapter"] == "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    assert bool(params["viral_reference"]) is viral
    assert params["covariates"] == []
    # Three tank units per arm, four residual df.
    assert manifest["deseq2_design"]["design"] == "~condition"
    assert manifest["deseq2_design"]["n_libraries_or_units"] == 6
    assert manifest["deseq2_design"]["residual_df"] == 4
    assert manifest["reference"]["genome_assembly"] == "GCF_963853765.1"


@pytest.mark.parametrize("comparison_id", list(RUN))
def test_results_have_standard_errors(comparison_id):
    with open(STUDY_DIR / f"{STUDY_ID}_{comparison_id}_dge_standardized.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) > 25000
    assert all(r["lfcSE"] not in ("", "NA") for r in rows if r["log2FoldChange"] not in ("", "NA"))


def test_every_fetched_run_matches_ena_read_count(samplesheet):
    with open(STUDY_DIR / "read_source.tsv") as fh:
        source = {r["run_accession"]: r for r in csv.DictReader(fh, delimiter="\t")}
    used = {r["run_accession"]: r for r in samplesheet if r["hpc"] in ("0", "24")}
    assert set(source) == set(used) and len(source) == 24
    for run, r in source.items():
        assert r["source_url"] == f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"
        assert r["reads_verified"] == used[run]["read_count"]
        assert len(r["sra_md5"]) == 32 and len(r["fastq_gz_md5"]) == 32


def test_dispersion_qc_is_committed_and_interpreted(study):
    with open(STUDY_DIR / "replicate_dispersion_qc.tsv") as fh:
        by_group = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["condition"] for r in by_group} == {"3_0", "3_24", "5_0", "5_24"}
    assert all(r["passes"] == "True" for r in by_group)
    with open(STUDY_DIR / "replicate_dispersion_qc_by_family.tsv") as fh:
        by_family = list(csv.DictReader(fh, delimiter="\t"))
    assert len(by_family) == 8
    # Within-family values fall below the single-oyster floor; the registry must
    # say why (pooled libraries) rather than leave it as a silent failure.
    assert any(r["passes"] == "False" for r in by_family)
    assert "not technical replicates" in study["limitations"]
