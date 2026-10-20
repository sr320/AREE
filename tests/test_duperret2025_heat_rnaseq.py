"""Curation tests for DUPERRET2025_HEAT_RNASEQ.

RNA-seq of Duperret et al. 2025 (Sci Total Environ): 30 C vs 23 C before
OsHV-1 challenge, deposited as two ENA projects by temperature (PRJEB86646
30 C, PRJEB86525 23 C) and registered as one study before any reanalysis.
The 23 C arm is the same set of libraries as DUPERRET2025_STARVE_RNASEQ's fed
arm. These tests pin the merge, the design, that shared control, and the
curation decisions that are easy to erode later: tank-level fitting, the
0 hpc temperature contrast as primary, and no duplicate 23 C infection
contrast.
"""
from __future__ import annotations

import csv
import json
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_ID = "DUPERRET2025_HEAT_RNASEQ"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
STARVE_DIR = DATA_DIR / "studies" / "DUPERRET2025_STARVE_RNASEQ"
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


def test_design_is_families_by_temperatures_by_seven_times(samplesheet, provenance):
    assert provenance["bioprojects"] == ["PRJEB86646", "PRJEB86525"]
    assert provenance["n_runs"] == len(samplesheet) == 84
    cells = Counter((r["family"], r["treatment_code"], r["hpc"]) for r in samplesheet)
    assert set(cells) == {(f, c, t) for f in ("F11N", "F14R") for c in ("3", "4") for t in TIMES}
    assert set(cells.values()) == {3}
    for r in samplesheet:
        assert r["temp_c"] == {"3": "23", "4": "30"}[r["treatment_code"]]
        assert float(r["temperature"]) == float(r["temp_c"])
        assert r["tank"] == f"{r['temp_c']}C_t{r['replicate']}"
        assert r["tank_hpc"] == f"{r['tank']}_{r['hpc']}hpc"


def test_control_arm_is_the_starvation_study_fed_arm(samplesheet, provenance):
    with open(STARVE_DIR / "samplesheet.tsv") as fh:
        fed = {r["run_accession"] for r in csv.DictReader(fh, delimiter="\t") if r["diet_code"] == "3"}
    control = {r["run_accession"] for r in samplesheet if r["treatment_code"] == "3"}
    assert control == fed and len(control) == 42
    assert set(provenance["shared_control_arm"]["runs"]) == control
    assert provenance["shared_control_arm"]["also_in_study"] == "DUPERRET2025_STARVE_RNASEQ"


def test_comparisons_fit_three_tank_units_per_arm(samplesheet, study):
    comps = {c["comparison_id"]: c for c in study["comparisons"]}
    assert set(comps) == {"heat30_vs_23c_0hpc", "heat30_oshv1_24hpc_vs_0hpc"}
    assert [cid for cid, c in comps.items() if c.get("meta_analysis_primary")] == ["heat30_vs_23c_0hpc"]
    primary = comps["heat30_vs_23c_0hpc"]
    assert primary["stressor_standardized"] == "temperature"
    assert primary["phenotype"] == "disease_resistance"
    assert primary["resilience_classification"] == "resilience"
    assert comps["heat30_oshv1_24hpc_vs_0hpc"]["resilience_classification"] == "disease"
    for (control, treatment), unit in [(("3_0", "4_0"), "tank"), (("4_0", "4_24"), "tank_hpc")]:
        rows = [r for r in samplesheet if r["condition"] in (control, treatment)]
        arms = {}
        for r in rows:
            arms.setdefault(r[unit], set()).add(r["condition"])
        assert all(len(a) == 1 for a in arms.values())
        assert Counter(next(iter(a)) for a in arms.values()) == {control: 3, treatment: 3}
        assert set(Counter(r[unit] for r in rows).values()) == {2}
    for c in comps.values():
        assert c["sample_size"] == 12 and c["biological_replicates"] == 3


def test_no_duplicate_23c_infection_contrast(study):
    # The 23 C (fed) 24 vs 0 hpc contrast already belongs to the starvation study.
    starve = yaml.safe_load((STUDIES_DIR / "DUPERRET2025_STARVE_RNASEQ.yaml").read_text())
    assert "fed_oshv1_24hpc_vs_0hpc" in {c["comparison_id"] for c in starve["comparisons"]}
    for c in study["comparisons"]:
        assert not (c["control_condition"].startswith("23 C") and "24 hpc" in c["treatment_condition"])


def test_accessions_and_shared_control_are_declared(study):
    accs = study["accessions"]
    assert accs["bioproject"] == "PRJEB86646"
    assert "PRJEB86646" in accs["ena"] and "PRJEB86525" in accs["ena"]
    assert accs["doi"] == "10.1016/j.scitotenv.2025.180026"
    assert "must not be pooled as independent studies" in study["limitations"]
    assert "tank_level_replication" in study["quality_flags"]


def test_config_is_single_end_and_fits_primary_on_tanks():
    config = (REPO_ROOT / "config" / f"{STUDY_ID}.config").read_text()
    assert "single_end = true" in config
    assert "replicate_unit  = 'tank'" in config


def test_fastq_manifest_covers_every_run_single_end(samplesheet):
    with open(STUDY_DIR / "fastq_manifest.tsv") as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["run_accession"] for r in manifest} == {r["run_accession"] for r in samplesheet}
    assert len(manifest) == len(samplesheet)
    assert all(len(r["md5"]) == 32 for r in manifest)


# ---- 2026-10-18 reanalysis outputs -------------------------------------------

RUN = {
    "heat30_vs_23c_0hpc": (["3_0"], ["4_0"], "tank", False),
    "heat30_oshv1_24hpc_vs_0hpc": (["4_0"], ["4_24"], "tank_hpc", True),
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
    used = {r["run_accession"]: r for r in samplesheet if r["condition"] in ("3_0", "4_0", "4_24")}
    assert set(source) == set(used) and len(source) == 18
    for run, r in source.items():
        assert r["source_url"] == f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"
        assert r["reads_verified"] == used[run]["read_count"]
        assert len(r["sra_md5"]) == 32 and len(r["fastq_gz_md5"]) == 32


def test_shared_control_reads_are_the_starvation_study_files():
    # Same .sra bytes as DUPERRET2025_STARVE_RNASEQ used for its fed 0 hpc arm.
    def sra_md5(path):
        with open(path) as fh:
            return {r["run_accession"]: r["sra_md5"] for r in csv.DictReader(fh, delimiter="\t")}
    ours, theirs = sra_md5(STUDY_DIR / "read_source.tsv"), sra_md5(STARVE_DIR / "read_source.tsv")
    shared = set(ours) & set(theirs)
    assert len(shared) == 6
    assert all(ours[r] == theirs[r] for r in shared)


def test_dispersion_qc_is_committed_and_interpreted(study):
    with open(STUDY_DIR / "replicate_dispersion_qc.tsv") as fh:
        by_group = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["condition"] for r in by_group} == {"3_0", "4_0", "4_24"}
    assert all(r["passes"] == "True" for r in by_group)
    with open(STUDY_DIR / "replicate_dispersion_qc_by_family.tsv") as fh:
        by_family = list(csv.DictReader(fh, delimiter="\t"))
    assert len(by_family) == 6
    # Pooled libraries fall below the single-oyster floor within a family; the
    # registry must say why rather than leave it as a silent failure.
    assert any(r["passes"] == "False" for r in by_family)
    assert "not technical replicates" in study["limitations"]


def test_infection_at_30c_took_in_f11n_only(study):
    with open(STUDY_DIR / f"{STUDY_ID}_heat30_oshv1_24hpc_vs_0hpc_viral_load.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    infected = {r["sample_id"] for r in rows if r["virus_detected"] == "True"}
    assert infected == {f"P_F11N_4_24_{i}" for i in (1, 2, 3)}
    assert not any(r["excluded_from_de"] == "True" for r in rows)
    assert "infection took in one family only" in study["limitations"]
