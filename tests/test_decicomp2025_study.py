"""Curation tests for DECICOMP2025_POMS_AGE.

AREE's first real metabolomics study (MetaboLights MTBLS13463), registered
from deposited ISA-Tab metadata before any analysis. These tests pin the
design (2 families x 3 ages x 6 whole animals), the finding that the deposited
MAFs hold no abundance values, and the curation decisions that are easy to
erode later: no processed data claimed, the linked EM-seq project kept out of
`accessions`, and the age contrasts kept separate against a 4-month reference.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "DECICOMP2025_POMS_AGE.yaml"
STUDY_DIR = DATA_DIR / "studies" / "DECICOMP2025_POMS_AGE"
FAMILIES = {"F14R", "H2D"}
AGES = {"C3": 4, "C7": 16, "C8": 28}


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def samplesheet():
    with open(STUDY_DIR / "samplesheet.tsv") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def provenance():
    return json.loads((STUDY_DIR / "metabolights_provenance.json").read_text())


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_design_is_two_families_by_three_ages_by_six(samplesheet, provenance):
    assert provenance["accession"] == "MTBLS13463"
    assert provenance["n_samples"] == len(samplesheet) == 36
    cells = Counter((r["family"], r["cohort"]) for r in samplesheet)
    assert set(cells) == {(f, c) for f in FAMILIES for c in AGES}
    assert set(cells.values()) == {6}
    for r in samplesheet:
        assert int(r["age_months"]) == AGES[r["cohort"]]
        assert r["condition"] == f"age{r['age_months']}mo"
        assert r["sample_id"] == f"{r['family']}_{r['cohort']}_{r['replicate']}"


def test_every_sample_has_one_raw_file_per_platform(samplesheet):
    for platform in ("GM", "NUCs", "TCA"):
        files = [r[f"raw_file_{platform}"] for r in samplesheet]
        assert all(f.startswith("FILES/RAW_FILES/") for f in files)
        assert len(set(files)) == 36


def test_deposited_mafs_hold_no_abundances(provenance, study):
    checks = provenance["maf_checks"]
    assert {p: c["n_features"] for p, c in checks.items()} == {"GM": 111, "NUCs": 21, "TCA": 13}
    for c in checks.values():
        assert c["n_abundance_values"] == 0
        assert c["n_sample_columns_in_samplesheet"] == 36
        assert c["n_database_identifiers"] == 0
    da = study["data_availability"]
    assert da["processed_data_available"] is False
    assert da["raw_data_available"] is True
    assert study["analysis_mode"] == "raw_reanalysis"
    assert all(c["results_file"] is None for c in study["comparisons"])


def test_linked_emseq_project_is_context_not_an_accession(study):
    accs = " ".join(str(v) for v in study["accessions"].values() if v)
    assert "PRJEB105019" not in accs
    assert study["accessions"]["other"] == "MTBLS13463"
    assert "PRJEB105019" in study["citation"]
    assert study["accessions"]["doi"] is None


def test_age_contrasts_are_separate_against_4_month_reference(study):
    comps = {c["comparison_id"]: c for c in study["comparisons"]}
    assert set(comps) == {"age28mo_vs_age4mo_basal", "age16mo_vs_age4mo_basal"}
    assert [cid for cid, c in comps.items() if c.get("meta_analysis_primary")] == ["age28mo_vs_age4mo_basal"]
    for c in comps.values():
        assert c["stressor_standardized"] == "pathogen_challenge"
        assert c["phenotype"] == "disease_resistance"
        assert c["phenotype_direction"] == "increase"
        assert c["resilience_classification"] == "resilience"
        assert c["tissue"] == "whole_animal"
        assert "4-month-old" in c["control_condition"]
        assert "not an unexposed control" in c["control_condition"]
        assert "inferred" in c["exposure_timing"]
        assert c["sample_size"] == 24 and c["biological_replicates"] == 6
    assert "ambiguous_phenotype_definition" in study["quality_flags"]


def test_samplesheet_checksum_matches_provenance(provenance):
    rel = "data/studies/DECICOMP2025_POMS_AGE/samplesheet.tsv"
    digest = hashlib.sha256((REPO_ROOT / rel).read_bytes()).hexdigest()
    assert provenance["outputs"][rel] == digest
