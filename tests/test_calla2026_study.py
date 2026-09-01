"""Tests for CALLA2026_OSHV — AREE's first registered raw_reanalysis study.

Nothing has been harmonized from this study yet; the reanalysis has not been
run. What these tests protect is the curation: that the registered design still
matches the design deposited in ENA, and that the study is not quietly
over-claimed as resilience evidence it does not contain.
"""
from __future__ import annotations

import csv
import json

import pytest
import yaml

from common import DATA_DIR, STUDIES_DIR, load_vocab

STUDY_PATH = STUDIES_DIR / "CALLA2026_OSHV.yaml"
STUDY_DIR = DATA_DIR / "studies" / "CALLA2026_OSHV"


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


# --------------------------------------------------------------------------- #
# The registered design must match the deposited design
# --------------------------------------------------------------------------- #


def test_study_validates(study):
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_samplesheet_matches_the_deposited_run_count(samplesheet, ena_provenance):
    assert len(samplesheet) == 42
    assert ena_provenance["n_runs"] == 42
    assert ena_provenance["bioproject"] == "PRJNA1329250"


def test_design_is_two_populations_by_four_viral_strains(samplesheet):
    breeds = {r["breed"] for r in samplesheet}
    strains = {r["viral_strain"] for r in samplesheet}
    assert breeds == {"Midori", "Myagi"}
    assert strains == {"Control", "Australia", "France", "USA"}


def test_every_group_is_replicated_enough_for_differential_expression(ena_provenance):
    """The screen that PRJNA623063 failed. n>=3 per group, or DESeq2 is meaningless."""
    assert ena_provenance["min_group_size"] >= 3
    assert len(ena_provenance["group_sizes"]) == 8


def test_each_comparison_matches_its_group_sizes(study, ena_provenance):
    """sample_size and biological_replicates must come from the deposited data."""
    groups = ena_provenance["group_sizes"]
    for comp in study["comparisons"]:
        population = "Midori" if comp["comparison_id"].startswith("midori") else "Myagi"
        strain = {"australia": "Australia", "france": "France", "usa": "USA"}[
            comp["comparison_id"].split("_")[2]
        ]
        treated = groups[f"{population}_{strain}"]
        control = groups[f"{population}_Control"]
        assert comp["biological_replicates"] == treated, comp["comparison_id"]
        assert comp["sample_size"] == treated + control, comp["comparison_id"]


def test_all_six_comparisons_are_registered(study):
    ids = {c["comparison_id"] for c in study["comparisons"]}
    assert len(ids) == 6
    assert all(("midori" in i) or ("miyagi" in i) for i in ids)


def test_fastq_manifest_covers_every_run_as_paired_end(samplesheet):
    with open(STUDY_DIR / "fastq_manifest.tsv") as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    runs_in_sheet = {r["run_accession"] for r in samplesheet}
    runs_in_manifest = {r["run_accession"] for r in manifest}
    assert runs_in_sheet == runs_in_manifest
    # Paired-end: two files per run, each with an ENA-supplied checksum.
    assert len(manifest) == 2 * len(runs_in_sheet)
    assert all(len(r["md5"]) == 32 for r in manifest)


# --------------------------------------------------------------------------- #
# The study must not be over-claimed
# --------------------------------------------------------------------------- #


def test_not_claimed_as_resilience_evidence(study):
    """The BioProject title says "tolerance", but no survival phenotype was
    measured for these animals. Registering this as disease_resistance would
    assert a resilience outcome the data does not contain."""
    for comp in study["comparisons"]:
        assert comp["phenotype"] == "disease_susceptibility", comp["comparison_id"]
        assert comp["resilience_classification"] == "disease", comp["comparison_id"]
        assert comp["phenotype_unit"] is None, (
            f"{comp['comparison_id']} claims a phenotype unit, but no phenotype "
            "was measured for these animals"
        )


def test_phenotype_classification_matches_the_ontology_default(study):
    """A mismatch would be a curator override needing written justification."""
    vocab = load_vocab("phenotype_ontology")
    for comp in study["comparisons"]:
        expected = vocab[comp["phenotype"]]["resilience_relevance"]
        assert comp["resilience_classification"] == expected


def test_unknown_exposure_parameters_are_null_not_invented(study):
    """Dose, duration, and timing are absent from the deposited metadata and the
    paper is paywalled. They must stay empty rather than be filled with a guess."""
    for comp in study["comparisons"]:
        assert comp["exposure_intensity"] is None
        assert comp["exposure_duration"] is None
        assert comp["exposure_timing"] is None


def test_analysis_status_reflects_partial_full_depth_reanalysis(study):
    """One full-depth comparison has been run; the remaining five are still pending."""
    assert study["analysis_status"] == "in_progress"
    assert study["qc_status"] == "in_progress"
    completed = {
        c["comparison_id"]: c["results_file"]
        for c in study["comparisons"]
        if c["results_file"] is not None
    }
    assert completed == {
        "miyagi_oshv1_usa_vs_control": (
            "data/studies/CALLA2026_OSHV/"
            "CALLA2026_OSHV_miyagi_oshv1_usa_vs_control_dge_standardized.tsv"
        )
    }


def test_ambiguous_phenotype_is_flagged(study):
    assert "ambiguous_phenotype_definition" in study["quality_flags"]


def test_reanalysis_target_is_the_crosswalk_annotation(study):
    """Reanalyzing against the annotation the crosswalk is built from is what
    spares this study the assembly crossing HESSER2024_VCOR carries."""
    from common import resolve_assembly

    assembly = resolve_assembly(study["genome_assembly"])
    assert assembly["assembly_id"] == "xbMagGiga1.1"
    assert assembly["is_ncbi_reference"] is True
    assert study["annotation_version"] == assembly["annotation_release"]


# --------------------------------------------------------------------------- #
# The BioProject cross-check
# --------------------------------------------------------------------------- #


def test_bioproject_mismatch_is_caught(tmp_path, monkeypatch):
    """A study YAML must not drift away from the data its sample sheet came from."""
    import intake.schema_validate as sv

    doc = yaml.safe_load(STUDY_PATH.read_text())
    doc["accessions"]["bioproject"] = "PRJNA9999999"
    path = tmp_path / "CALLA2026_OSHV.yaml"
    path.write_text(yaml.safe_dump(doc))
    monkeypatch.setattr(sv, "DATA_DIR", DATA_DIR)

    result = sv.validate_study_file(path)
    assert not result.valid
    assert any("PRJNA9999999" in e for e in result.errors)
