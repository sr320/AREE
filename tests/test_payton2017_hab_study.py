"""Curation tests for PAYTON2017_HAB_AMINUTUM.

Payton et al. 2017 (PRJEB18614 / E-MTAB-5345): adult gill after 6 days on
the toxic dinoflagellate Alexandrium minutum vs non-toxic Heterocapsa
triquetra, 13 pooled libraries per diet (one per sampling time), each on two
lanes. AREE's first harmful-algal-bloom study. These tests pin the diet
labels against the ArrayExpress SDRF, the lane structure (two runs per
library, summed by replicate_unit), and the classification.
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "PAYTON2017_HAB_AMINUTUM.yaml"
STUDY_DIR = DATA_DIR / "studies" / "PAYTON2017_HAB_AMINUTUM"
CONFIG = REPO_ROOT / "config" / "PAYTON2017_HAB_AMINUTUM.config"


def _tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def samplesheet():
    return _tsv(STUDY_DIR / "samplesheet.tsv")


@pytest.fixture(scope="module")
def sdrf_by_run():
    rows = _tsv(STUDY_DIR / "E-MTAB-5345.sdrf.txt")
    return {r["Comment[ENA_RUN]"]: r for r in rows}


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_diet_matches_the_arrayexpress_sdrf_for_every_run(samplesheet, sdrf_by_run):
    assert set(sdrf_by_run) == {r["run_accession"] for r in samplesheet}
    for row in samplesheet:
        assert row["growth_condition"] == sdrf_by_run[row["run_accession"]]["Factor Value[feeding condition]"]


def test_two_lanes_per_library_thirteen_libraries_per_diet(samplesheet, sdrf_by_run):
    runs_per_library = Counter(r["sample_accession"] for r in samplesheet)
    assert set(runs_per_library.values()) == {2}
    libraries = defaultdict(set)
    times = defaultdict(set)
    for row in samplesheet:
        libraries[row["condition"]].add(row["sample_accession"])
        times[row["condition"]].add(sdrf_by_run[row["run_accession"]]["Factor Value[time]"])
    assert {k: len(v) for k, v in libraries.items()} == {"Alexandrium_minutum": 13, "Heterocapsa_triquetra": 13}
    # Time-matched: both diets sampled at the same 13 times.
    assert times["Alexandrium_minutum"] == times["Heterocapsa_triquetra"]
    assert len(times["Alexandrium_minutum"]) == 13


def test_lanes_are_summed_so_libraries_are_the_replicates(study):
    assert re.search(r"replicate_unit\s*=\s*'sample_accession'", CONFIG.read_text())
    (comp,) = study["comparisons"]
    assert comp["biological_replicates"] == 13 and comp["sample_size"] == 26


def test_hab_exposure_is_stress_response_with_its_replication_caveat(study):
    (comp,) = study["comparisons"]
    assert comp["stressor_standardized"] == "harmful_algal_bloom"
    assert comp["phenotype"] == "toxicant_tolerance"
    assert comp["resilience_classification"] == "stress_response"
    assert comp["meta_analysis_primary"] is True
    assert "tank_level_replication" in study["quality_flags"]


def test_provenance_and_manifest(samplesheet):
    prov = json.loads((STUDY_DIR / "ena_provenance.json").read_text())
    assert prov["bioproject"] == "PRJEB18614" and prov["n_runs"] == 52
    manifest = _tsv(STUDY_DIR / "fastq_manifest.tsv")
    assert len(manifest) == 104
    assert all(len(r["md5"]) == 32 for r in manifest)
