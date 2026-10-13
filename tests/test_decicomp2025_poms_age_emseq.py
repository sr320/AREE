"""Curation tests for DECICOMP2025_POMS_AGE_EMSEQ.

EM-seq arm of Ifremer's DECICOMP POMS age experiment (ENA PRJEB105019),
registered from deposited metadata before any reanalysis. These tests pin the
design parsed from experiment_title and its agreement with the BioSample
attributes, and the curation decisions that are easy to erode later:
individual oysters as replicates, the basal age contrast as primary (pairing
with DECICOMP2025_POMS_AGE), the R/S contrast flagged as ambiguous, and the
metabolomics arm kept out of `accessions`.
"""
from __future__ import annotations

import csv
import json
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_ID = "DECICOMP2025_POMS_AGE_EMSEQ"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
GROUP_PHENOTYPE = {"T0": "No infected (Basal)", "R": "Resistant", "S": "Susceptible"}
COHORT_AGE = {"C3": "4 months", "C7": "16 months", "C8": "28 months"}


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


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_design_is_five_groups_of_six_per_family(samplesheet, ena_provenance):
    assert ena_provenance["bioproject"] == "PRJEB105019"
    assert ena_provenance["n_runs"] == len(samplesheet) == 60
    assert ena_provenance["library_layouts"] == ["PAIRED"]
    cells = Counter((r["condition"], r["family"]) for r in samplesheet)
    groups = {"C3_T0", "C7_T0", "C8_T0", "C8_R", "C8_S"}
    assert set(cells) == {(g, f) for g in groups for f in ("F14R", "H2D")}
    assert set(cells.values()) == {6}


def test_parsed_labels_agree_with_biosample_attributes(samplesheet):
    for r in samplesheet:
        assert r["condition"] == f"{r['cohort']}_{r['group']}"
        assert r["phenotype"] == GROUP_PHENOTYPE[r["group"]]
        assert r["age"] == COHORT_AGE[r["cohort"]]
        assert r["isolate"] == f"Family {r['family']}"
        assert r["sample_id"] == f"{r['family']}_{r['cohort']}_{r['group']}_S{r['replicate']}"


def test_methylkit_sheet_codes_primary_contrast_by_run(samplesheet):
    with open(STUDY_DIR / "methylkit_samplesheet.csv") as fh:
        coded = {r["sample_id"]: r["treatment"] for r in csv.DictReader(fh)}
    expected = {r["run_accession"]: "1" if r["condition"] == "C8_T0" else "0"
                for r in samplesheet if r["condition"] in ("C3_T0", "C8_T0")}
    assert coded == expected


def test_comparisons_and_primary(study):
    comps = {c["comparison_id"]: c for c in study["comparisons"]}
    assert set(comps) == {"age28mo_vs_age4mo_basal", "age16mo_vs_age4mo_basal",
                          "resistant_vs_susceptible_28mo_500hpc"}
    assert [cid for cid, c in comps.items() if c.get("meta_analysis_primary")] == ["age28mo_vs_age4mo_basal"]
    for c in comps.values():
        assert c["stressor_standardized"] == "pathogen_challenge"
        assert c["phenotype"] == "disease_resistance"
        assert c["resilience_classification"] == "resilience"
        assert c["sample_size"] == 24 and c["biological_replicates"] == 6
        assert c["results_file"] is None
    assert "criteria not deposited" in comps["resistant_vs_susceptible_28mo_500hpc"]["phenotype_unit"]
    assert "ambiguous_phenotype_definition" in study["quality_flags"]


def test_individuals_are_replicates_no_replicate_unit():
    config = (REPO_ROOT / "config" / f"{STUDY_ID}.config").read_text()
    params = config.split("params {", 1)[1]
    assert "replicate_unit" not in params.split("process {", 1)[0]
    assert "design_sheet" not in params


def test_metabolomics_arm_is_context_not_an_accession(study):
    accs = " ".join(str(v) for v in study["accessions"].values() if v)
    assert "MTBLS13463" not in accs
    assert "MTBLS13463" in study["citation"]
    assert study["accessions"]["bioproject"] == "PRJEB105019"
    assert study["accessions"]["doi"] is None


def test_fastq_manifests_cover_every_run_paired_with_md5s(samplesheet):
    rows = []
    for name in ("fastq_manifest.tsv", "fastq_manifest_supplement.tsv"):
        with open(STUDY_DIR / name) as fh:
            rows += [(name, r) for r in csv.DictReader(fh, delimiter="\t")]
    per_run = Counter(r["run_accession"] for _, r in rows)
    assert set(per_run) == {r["run_accession"] for r in samplesheet}
    assert set(per_run.values()) == {2}
    assert all(len(r["md5"]) == 32 and int(r["bytes"]) > 0 for _, r in rows)
    # ERR15996548 has no ENA-processed FASTQ; its submitted files are the only source.
    supplement = {r["run_accession"] for name, r in rows if name == "fastq_manifest_supplement.tsv"}
    assert supplement == {"ERR15996548"}
