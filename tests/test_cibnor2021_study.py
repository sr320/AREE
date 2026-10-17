"""Curation tests for CIBNOR2021_HEAT_RRSS.

AREE's first real methylation study, and its first real thermal-resilience
contrast. The project deposits all 12 libraries under one BioSample, so the
design is parsed from experiment_title. These tests pin that parse, the family
structure behind it, and the decisions that are easy to erode later:
biological replicates counted as families, no borrowed DOI, SS as the
reference, and a result tested on 2 vs 2 families that is reported as the
null it is.
"""
from __future__ import annotations

import csv
import json

import pytest
import yaml

from common import DATA_DIR, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "CIBNOR2021_HEAT_RRSS.yaml"
STUDY_DIR = DATA_DIR / "studies" / "CIBNOR2021_HEAT_RRSS"
FAMILIES = {"RR": {"52", "59"}, "SS": {"05", "35"}}


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


def test_design_is_two_families_per_class_three_libraries_each(samplesheet, ena_provenance):
    assert ena_provenance["bioproject"] == "PRJNA690951"
    assert ena_provenance["library_strategies"] == ["Bisulfite-Seq"]
    assert ena_provenance["n_runs"] == len(samplesheet) == 12
    # One BioSample for every library: the reason run labels were needed.
    assert ena_provenance["n_samples"] == 1
    assert ena_provenance["group_sizes"] == {"RR": 6, "SS": 6}

    for cls, fams in FAMILIES.items():
        rows = [r for r in samplesheet if r["condition"] == cls]
        assert {r["family"] for r in rows} == fams
        for fam in fams:
            assert sorted(r["replicate"] for r in rows if r["family"] == fam) == ["1", "2", "3"]


def test_parsed_labels_are_unique_and_recorded(samplesheet, ena_provenance):
    labels = ena_provenance["run_labels"]["labels"]
    assert ena_provenance["run_labels"]["field"] == "experiment_title"
    assert len(set(labels.values())) == 12
    for row in samplesheet:
        assert labels[row["run_accession"]] == row["sample_id"]
        assert row["sample_id"] == f"{row['condition']}{row['family']}-R{row['replicate']}"


def test_methylkit_sheet_codes_ss_as_reference_by_run_accession(samplesheet):
    with open(STUDY_DIR / "methylkit_samplesheet.csv") as fh:
        coded = {r["sample_id"]: r["treatment"] for r in csv.DictReader(fh)}
    expected = {r["run_accession"]: "1" if r["condition"] == "RR" else "0" for r in samplesheet}
    assert coded == expected


def test_replication_is_counted_in_families_not_libraries(study):
    (comp,) = study["comparisons"]
    assert comp["sample_size"] == 12
    assert comp["biological_replicates"] == 2
    assert comp["meta_analysis_primary"] is True
    for flag in ("family_level_replication", "low_replication", "ambiguous_phenotype_definition"):
        assert flag in study["quality_flags"]


def test_tolerance_contrast_under_one_challenge_not_an_exposure(study):
    (comp,) = study["comparisons"]
    assert comp["stressor_standardized"] == "temperature"
    assert comp["phenotype"] == "thermal_tolerance"
    assert comp["resilience_classification"] == "resilience"
    assert comp["phenotype_direction"] == "increase"
    assert comp["tissue"] == "gill"
    assert "not an unexposed control" in comp["control_condition"]


def test_rnaseq_paper_is_context_not_the_source_doi(study):
    assert study["accessions"]["doi"] is None
    assert "10.1016/j.cbd.2023.101089" in study["citation"]
    assert "does not analyze these data" in study["citation"]


def test_fastq_manifest_covers_every_run_as_paired_end(samplesheet):
    with open(STUDY_DIR / "fastq_manifest.tsv") as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    runs_in_sheet = {r["run_accession"] for r in samplesheet}
    assert {r["run_accession"] for r in manifest} == runs_in_sheet
    assert len(manifest) == 2 * len(runs_in_sheet)
    assert all(len(r["md5"]) == 32 for r in manifest)



def test_config_tests_families_with_limma_on_pooled_tiles():
    from common import REPO_ROOT

    config = (REPO_ROOT / "config" / "CIBNOR2021_HEAT_RRSS.config").read_text()
    for setting in ("replicate_unit = 'family'", "dmr_engine        = 'limma'",
                    "bismark_align_args = '--score_min L,0,-0.6'", "dmr_features      = 'genes'",
                    "extractor_args     = '--ignore_r2 20 --ignore_3prime 1 --ignore_3prime_r2 1'"):
        assert setting in config, setting


def test_manifest_records_a_family_level_limma_run():
    stem = "CIBNOR2021_HEAT_RRSS_rr_vs_ss_day30_oscillating_heat"
    manifest = json.loads((STUDY_DIR / f"{stem}_workflow_manifest.json").read_text())
    params = manifest["parameters"]
    assert params["replicate_unit"] == "family"
    assert params["dmr_engine"] == "limma"
    assert params["unit_min_coverage"] == 10
    with open(STUDY_DIR / f"{stem}_dmr_qc.tsv") as fh:
        (qc,) = list(csv.DictReader(fh, delimiter="\t"))
    assert qc["n_replicate_units"] == "4"
    assert qc["engine"] == "limma"


def test_result_is_a_low_power_null_not_a_list_of_dmrs(study):
    """Gene-body units plus intergenic tiles, none at q < 0.05: a null, with effect sizes."""
    import pandas as pd

    from common import REPO_ROOT

    (comp,) = study["comparisons"]
    df = pd.read_csv(REPO_ROOT / comp["results_file"], sep="\t")
    assert len(df) == 132587
    genes = df[df["region_id"].str.startswith("GENE_")]
    assert len(genes) == 31166 and genes["gene_id"].is_unique  # one record per gene
    assert (df["qvalue"] < 0.05).sum() == 0
    # Poolable: every tile has a standard error and an unadjusted p-value.
    assert df["se_meth_diff_percent"].notna().all() and df["pvalue"].notna().all()
    assert set(df["analysis_method"]) == {"limma_pooled_genes"}
    assert study["qc_status"] == "passed_with_warnings"
    assert "null at very low power" in STUDY_PATH.read_text()


def test_libraries_vary_like_distinct_animals_within_each_family():
    with open(STUDY_DIR / "replicate_dispersion_qc.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["family"] for r in rows} == FAMILIES["RR"] | FAMILIES["SS"]
    # Technical replicates would sit near 1 (binomial); these are 6.7-8.1.
    assert all(float(r["median_between_library_variance_over_binomial"]) > 3 for r in rows)


def test_every_library_mapped_and_converted_well():
    with open(STUDY_DIR / "library_qc.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 12
    assert all(float(r["unique_mapping_pct"]) >= 50 for r in rows)
    assert all(float(r["chh_methylation_pct"]) <= 1.5 for r in rows)
