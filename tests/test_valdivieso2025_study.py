"""Curation tests for VALDIVIESO2025_POMS_RS_EMSEQ.

EM-seq of gill and mantle biopsies taken from the same 10 oysters before (T0)
and after (T1) an experimental POMS challenge, each oyster classed resistant or
susceptible by its own outcome (ENA PRJEB81880; Valdivieso et al. 2025).
Registered before any reanalysis. These tests pin the deposited design and the
curation decisions that are easy to erode later: individuals (not families) as
replicates, one resistant and one susceptible animal per population, the
contrast type of each comparison, and gill as the primary within each type.
"""
from __future__ import annotations

import csv
import json
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_ID = "VALDIVIESO2025_POMS_RS_EMSEQ"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
PHENOTYPE = {"R": "resistant", "S": "susceptible"}
TISSUE = {"G": "gills", "M": "mantle"}
POPULATIONS = {f"GT{i}" for i in range(1, 6)}

EXPECTED_CONTRASTS = {
    "resistant_vs_susceptible_gill_t0": ("baseline_phenotype_difference", True),
    "resistant_vs_susceptible_mantle_t0": ("baseline_phenotype_difference", False),
    "resistant_vs_susceptible_gill_t1": ("phenotype_difference_under_exposure", True),
    "resistant_vs_susceptible_mantle_t1": ("phenotype_difference_under_exposure", False),
    "susceptible_t1_vs_t0_gill": ("exposure_response", True),
    "resistant_t1_vs_t0_gill": ("exposure_response", False),
    "susceptible_t1_vs_t0_mantle": ("exposure_response", False),
    "resistant_t1_vs_t0_mantle": ("exposure_response", False),
}


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


def test_design_is_eight_cells_of_five(samplesheet, ena_provenance):
    assert ena_provenance["bioproject"] == "PRJEB81880"
    assert ena_provenance["n_runs"] == len(samplesheet) == 40
    assert ena_provenance["library_layouts"] == ["PAIRED"]
    assert ena_provenance["library_strategies"] == ["Bisulfite-Seq"]
    cells = Counter(r["condition"] for r in samplesheet)
    assert set(cells) == {f"{t}_{p}_{tp}" for t in "GM" for p in "RS" for tp in ("T0", "T1")}
    assert set(cells.values()) == {5}
    # Each cell holds one library from each of the five populations.
    for cond in cells:
        assert {r["population"] for r in samplesheet if r["condition"] == cond} == POPULATIONS


def test_ten_animals_one_resistant_one_susceptible_per_population(samplesheet):
    animals = Counter((r["population"], r["code"], r["phenotype_class"]) for r in samplesheet)
    assert len(animals) == 10
    # Each animal: gill and mantle at T0 and T1.
    assert set(animals.values()) == {4}
    for key in animals:
        libs = {(r["tissue_code"], r["time"]) for r in samplesheet
                if (r["population"], r["code"], r["phenotype_class"]) == key}
        assert libs == {(t, tp) for t in "GM" for tp in ("T0", "T1")}
    per_pop = Counter((pop, cls) for pop, _, cls in animals)
    assert set(per_pop) == {(pop, cls) for pop in POPULATIONS for cls in "RS"}
    assert set(per_pop.values()) == {1}


def test_parsed_labels_agree_with_biosample_attributes(samplesheet):
    for r in samplesheet:
        assert r["condition"] == f"{r['tissue_code']}_{r['phenotype_class']}_{r['time']}"
        assert r["phenotype"] == PHENOTYPE[r["phenotype_class"]]
        assert r["tissue_type"] == TISSUE[r["tissue_code"]]
        assert r["age"] == "10 months old"
        assert f"Family {r['population'][2]} " in r["sample_description"]
        assert f"({r['time']})" in r["sample_description"]
        assert f"({r['code']})" in r["sample_description"]
        assert r["sample_id"] == (
            f"{r['population']}_{r['code']}_{r['time']}_{r['tissue_code']}_{r['phenotype_class']}"
        )


def test_methylkit_sheet_codes_primary_baseline_by_run(samplesheet):
    with open(STUDY_DIR / "methylkit_samplesheet.csv") as fh:
        coded = {r["sample_id"]: r["treatment"] for r in csv.DictReader(fh)}
    expected = {r["run_accession"]: "1" if r["condition"] == "G_R_T0" else "0"
                for r in samplesheet if r["condition"] in ("G_S_T0", "G_R_T0")}
    assert coded == expected


def test_comparisons_contrast_types_and_primaries(study):
    comps = {c["comparison_id"]: c for c in study["comparisons"]}
    assert set(comps) == set(EXPECTED_CONTRASTS)
    for cid, (contrast, primary) in EXPECTED_CONTRASTS.items():
        c = comps[cid]
        assert c["contrast_type"] == contrast, cid
        assert bool(c.get("meta_analysis_primary")) is primary, cid
        assert c["stressor_standardized"] == "pathogen_challenge"
        # Individuals are the replicates: 5 per arm, 10 libraries per comparison.
        assert c["sample_size"] == 10 and c["biological_replicates"] == 5
        assert c["results_file"] is None
    # Exactly one primary per contrast type, and it is the gill comparison.
    primaries = Counter(c["contrast_type"] for c in comps.values() if c.get("meta_analysis_primary"))
    assert primaries == Counter({t: 1 for t, _ in EXPECTED_CONTRASTS.values()})
    assert all(comps[cid]["tissue"] == "gill"
               for cid, (_, primary) in EXPECTED_CONTRASTS.items() if primary)


def test_phenotype_and_classification_by_contrast(study):
    for c in study["comparisons"]:
        if c["contrast_type"] == "exposure_response":
            assert c["phenotype"] == "disease_susceptibility"
            assert c["resilience_classification"] == "disease"
        else:
            # Resistant is the treatment arm, susceptible the reference.
            assert c["phenotype"] == "disease_resistance"
            assert c["resilience_classification"] == "resilience"
            assert c["phenotype_direction"] == "increase"
            assert "resistant" in c["treatment_condition"]
            assert "susceptible" in c["control_condition"]


def test_registered_but_not_run(study):
    assert study["analysis_mode"] == "raw_reanalysis"
    assert study["analysis_status"] == "not_started"
    assert study["qc_status"] == "not_started"
    assert study["accessions"]["doi"] == "10.1016/j.scitotenv.2025.178385"
    # The phenotype is individual, not family-level.
    assert "family_level_replication" not in study["quality_flags"]


def test_individuals_are_replicates_no_replicate_unit():
    config = (REPO_ROOT / "config" / f"{STUDY_ID}.config").read_text()
    params = config.split("params {", 1)[1].split("process {", 1)[0]
    assert "replicate_unit" not in params
    assert "design_sheet" not in params
    assert "resistant_vs_susceptible_gill_t0" in params


def test_fastq_manifest_covers_every_run_paired_with_md5s(samplesheet):
    with open(STUDY_DIR / "fastq_manifest.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    per_run = Counter(r["run_accession"] for r in rows)
    assert set(per_run) == {r["run_accession"] for r in samplesheet}
    assert set(per_run.values()) == {2}
    assert all(len(r["md5"]) == 32 and int(r["bytes"]) > 0 for r in rows)
    assert {r["source"] for r in rows} == {"fastq_ftp"}
