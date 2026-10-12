"""Curation tests for LI2022_HEAT_PLOIDY.

Diploid and triploid gill at time 0, after 12 h of acute heat, and at the end
of prolonged heat at 28 C, n=3 per group, reanalyzed from raw reads. These
tests pin the design read from ENA and the decisions that are easy to erode
later: one animal per library, the four prespecified heat-vs-baseline
comparisons with diploid prolonged heat as the only pooled one, exposure_only
classification, exposure details left null rather than guessed, the
publication link recorded as inferred, the config fix that stops every
comparison being labelled as the primary, and the borderline t_f dispersion.
"""
from __future__ import annotations

import csv
import json
import re

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_ID = "LI2022_HEAT_PLOIDY"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
PRIMARY = "diploid_prolonged_heat_vs_baseline"
# comparison_id -> (control level, treatment level), prespecified 2026-10-10.
COMPARISONS = {
    "diploid_prolonged_heat_vs_baseline": ("d_0", "d_f"),
    "diploid_acute_heat_12h_vs_baseline": ("d_0", "d_12"),
    "triploid_prolonged_heat_vs_baseline": ("t_0", "t_f"),
    "triploid_acute_heat_12h_vs_baseline": ("t_0", "t_12"),
}
# Genes at padj < 0.05 in each committed result (Pori, 2026-10-10/11).
DEG_COUNTS = {
    "diploid_prolonged_heat_vs_baseline": 849,
    "diploid_acute_heat_12h_vs_baseline": 2522,
    "triploid_prolonged_heat_vs_baseline": 2209,
    "triploid_acute_heat_12h_vs_baseline": 971,
}


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def comparisons(study):
    return {c["comparison_id"]: c for c in study["comparisons"]}


@pytest.fixture(scope="module")
def samplesheet():
    with open(STUDY_DIR / "samplesheet.tsv") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def ena_provenance():
    return json.loads((STUDY_DIR / "ena_provenance.json").read_text())


def _manifest(comparison_id: str) -> dict:
    path = STUDY_DIR / f"{STUDY_ID}_{comparison_id}_workflow_manifest.json"
    return json.loads(path.read_text())


# --------------------------------------------------------------------------- #
# The registered design must match the deposited design
# --------------------------------------------------------------------------- #


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_samplesheet_matches_the_deposit(samplesheet, ena_provenance):
    assert ena_provenance["bioproject"] == "PRJNA762441"
    assert ena_provenance["n_runs"] == ena_provenance["n_samples"] == len(samplesheet) == 18
    assert ena_provenance["library_strategies"] == ["RNA-Seq"]
    assert ena_provenance["library_layouts"] == ["PAIRED"]
    assert ena_provenance["group_sizes"] == dict.fromkeys(("d_0", "d_12", "d_f", "t_0", "t_12", "t_f"), 3)
    assert {r["tissue"] for r in samplesheet} == {"Gills"}


def test_condition_is_parsed_from_the_sample_alias(samplesheet):
    """sample_alias A<time><ploidy>_<rep> is the only place the design lives."""
    for row in samplesheet:
        time, ploidy, rep = re.fullmatch(r"A(0|12|f)(d|t)_([1-3])", row["sample_alias"]).groups()
        assert row["condition"] == f"{ploidy}_{time}"
        assert (row["ploidy"], row["timepoint"], row["replicate"]) == (ploidy, time, rep)


def test_every_library_is_a_different_oyster(samplesheet):
    """Distinct BioSamples, and isolation_source numbers oysters 1-9 per ploidy,
    three per time point. Lethal gill sampling means one animal per library."""
    assert len({r["sample_accession"] for r in samplesheet}) == 18
    oyster_numbers = {"0": {1, 2, 3}, "12": {4, 5, 6}, "f": {7, 8, 9}}
    seen: dict[str, set[int]] = {"d": set(), "t": set()}
    for row in samplesheet:
        word, number = re.fullmatch(r"(diploid|triploid) Crassostrea gigas([1-9])", row["isolation_source"]).groups()
        assert word[0] == row["ploidy"]
        assert int(number) in oyster_numbers[row["timepoint"]]
        seen[row["ploidy"]].add(int(number))
    assert seen == {"d": set(range(1, 10)), "t": set(range(1, 10))}


# --------------------------------------------------------------------------- #
# Prespecified comparisons and how they are classified
# --------------------------------------------------------------------------- #


def test_four_prespecified_heat_vs_baseline_comparisons(comparisons):
    assert set(comparisons) == set(COMPARISONS)
    for cid, (control, treatment) in COMPARISONS.items():
        comp = comparisons[cid]
        # Each arm names its sample-sheet level, e.g. "(d_0; A0d_1-3)".
        assert re.search(rf"[ (]{control};", comp["control_condition"])
        assert re.search(rf"[ (]{treatment};", comp["treatment_condition"])
        # Each arm stays within one ploidy: heat vs that ploidy's own baseline.
        assert control[0] == treatment[0] == cid[0]
        assert control.endswith("_0")
        assert comp["sample_size"] == 6
        assert comp["biological_replicates"] == 3


def test_only_diploid_prolonged_heat_is_pooled(comparisons):
    primary = [cid for cid, c in comparisons.items() if c.get("meta_analysis_primary")]
    assert primary == [PRIMARY]


def test_every_comparison_is_exposure_only(comparisons):
    """Survival was compared between ploidy groups, not measured on these 18
    animals. Do not upgrade to `resilience` on the strength of the paper."""
    for comp in comparisons.values():
        assert comp["stressor_standardized"] == "temperature"
        assert comp["phenotype"] == "thermal_tolerance"
        assert comp["resilience_classification"] == "exposure_only"
        assert comp["phenotype_direction"] == "not_applicable"
        assert comp["tissue"] == "gill"


def test_unread_exposure_details_are_null_not_guessed(comparisons):
    """Until the methods are read: no prolonged duration, no acute temperature."""
    for cid, comp in comparisons.items():
        if "prolonged" in cid:
            assert comp["exposure_duration"] is None
            assert comp["exposure_intensity"].startswith("28 C")
        else:
            assert comp["exposure_intensity"] is None


def test_publication_link_is_recorded_as_inferred(study):
    assert study["accessions"]["doi"] == "10.1016/j.aquaculture.2022.738219"
    assert "INFERRED" in STUDY_PATH.read_text()
    assert "Confirm" in STUDY_PATH.read_text()


# --------------------------------------------------------------------------- #
# Run configuration and the committed results
# --------------------------------------------------------------------------- #


def test_config_lets_comparison_id_be_overridden():
    """params.rnaseq.comparison_id is read before --comparison_id. If it is set,
    every comparison is labelled as the primary (2026-10-11)."""
    config = (REPO_ROOT / "config" / f"{STUDY_ID}.config").read_text()
    rnaseq_block = re.search(r"rnaseq\s*\{(.*?)\}", config, re.S).group(1)
    assert re.search(r"^\s*comparison_id\s*=\s*null\s*$", rnaseq_block, re.M)
    params = dict(re.findall(r"^    (comparison_id|control_level|treatment_level)\s*=\s*'([^']*)'", config, re.M))
    assert params == {"comparison_id": PRIMARY, "control_level": "d_0", "treatment_level": "d_f"}


@pytest.mark.parametrize("comparison_id", sorted(COMPARISONS))
def test_manifest_matches_the_comparison_it_is_named_for(comparison_id):
    manifest = _manifest(comparison_id)
    control, treatment = COMPARISONS[comparison_id]
    assert manifest["study_id"] == STUDY_ID
    assert manifest["comparison_id"] == comparison_id
    assert manifest["mode"] == "raw_reanalysis"
    assert manifest["parameters"]["control_levels"] == [control]
    assert manifest["parameters"]["treatment_levels"] == [treatment]
    assert manifest["reference"]["genome_assembly"] == "GCF_963853765.1"


@pytest.mark.parametrize("comparison_id", sorted(COMPARISONS))
def test_standardized_result_is_poolable(comparisons, comparison_id):
    import pandas as pd

    path = comparisons[comparison_id]["results_file"]
    assert path == f"data/studies/{STUDY_ID}/{STUDY_ID}_{comparison_id}_dge_standardized.tsv"
    df = pd.read_csv(REPO_ROOT / path, sep="\t")
    assert {"gene_id", "log2FoldChange", "lfcSE", "pvalue", "padj"} <= set(df.columns)
    assert df["lfcSE"].notna().all()
    assert (df["padj"] < 0.05).sum() == DEG_COUNTS[comparison_id]


def test_replicates_are_biological_and_t_f_is_flagged(study):
    """All six groups pass; t_f (0.0095) is a borderline pass and must stay
    documented as one, not silently treated like the others."""
    with open(STUDY_DIR / "replicate_dispersion_qc.tsv") as fh:
        rows = {r["condition"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert set(rows) == {"d_0", "d_12", "d_f", "t_0", "t_12", "t_f"}
    assert all(r["passes"] == "True" for r in rows.values())
    assert all(float(r["median_dispersion"]) > 0.005 for r in rows.values())
    assert float(rows["t_f"]["median_dispersion"]) < 0.01
    assert "borderline" in study["provenance"]["curation_notes"]
    assert study["qc_status"] == "passed_with_warnings"
    assert study["analysis_status"] == "harmonized"
