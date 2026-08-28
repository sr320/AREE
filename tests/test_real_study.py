"""Tests covering behaviours that only appear once a real published study is registered.

Each of these encodes something that the synthetic demo data could not have
exercised, and most of them are regressions for bugs that the first real study
actually exposed.
"""
import json
from pathlib import Path

import pandas as pd
import pytest

from common import REPO_ROOT
from harmonize.identifiers import (
    DEMO_CROSSWALK_PATH,
    crosswalk_for_study,
    identifier_candidates,
    resolve_identifier,
    set_crosswalk_path,
)
from intake.schema_validate import validate_study_file

STUDY_ID = "HESSER2024_VCOR"
STUDY_YAML = REPO_ROOT / "registry" / "studies" / f"{STUDY_ID}.yaml"
STUDY_DIR = REPO_ROOT / "data" / "studies" / STUDY_ID
REAL_CROSSWALK = REPO_ROOT / "data" / "reference" / "crosswalk" / "mgigas_gene_id_crosswalk.tsv"


@pytest.fixture
def real_crosswalk():
    if not REAL_CROSSWALK.exists():
        pytest.skip("real crosswalk not built; run `aree build-crosswalk`")
    set_crosswalk_path(REAL_CROSSWALK)
    yield REAL_CROSSWALK
    set_crosswalk_path(None)


# --------------------------------------------------------------------------- #
# Registration
# --------------------------------------------------------------------------- #


def test_real_study_registration_is_valid():
    result = validate_study_file(STUDY_YAML)
    assert result.valid, result.errors


def test_real_study_is_not_marked_simulated():
    import yaml

    doc = yaml.safe_load(STUDY_YAML.read_text())
    assert doc["simulated"] is False
    assert doc["accessions"]["doi"], "a real study must carry at least one accession"


def test_real_study_records_its_limitations():
    """A significance-filtered source table without standard errors must say so."""
    import yaml

    doc = yaml.safe_load(STUDY_YAML.read_text())
    limitations = doc["limitations"].lower()
    assert "standard error" in limitations
    assert "processed_only" in doc["quality_flags"]
    assert "low_replication" in doc["quality_flags"]


def test_intake_provenance_matches_the_generated_result_files():
    prov = json.loads((STUDY_DIR / "intake_provenance.json").read_text())
    from common import sha256_file

    assert prov["source"]["license"].startswith("CC BY")
    for conv in prov["conversions"]:
        path = REPO_ROOT / conv["output_file"]
        assert path.exists(), path
        assert conv["output_sha256"] == sha256_file(path), (
            f"{path} changed since intake provenance was written; re-run the converter"
        )


# --------------------------------------------------------------------------- #
# Identifier decoration (the source reports 'gene-LOC123|LOC123')
# --------------------------------------------------------------------------- #


def test_identifier_candidates_strips_gff_decoration():
    cands = identifier_candidates("gene-LOC105320749|LOC105320749")
    assert cands[0] == "gene-LOC105320749|LOC105320749", "verbatim value must be tried first"
    assert "LOC105320749" in cands


def test_identifier_candidates_leaves_a_plain_id_alone():
    assert identifier_candidates("105320749") == ["105320749"]


def test_decorated_identifier_resolves_exactly(real_crosswalk):
    resolved = resolve_identifier("gene-LOC105320749|LOC105320749", "ncbi_gene_id")
    assert resolved.mapping_confidence == "exact"
    assert resolved.feature_id_standardized == "105320749"


def test_loc_form_and_numeric_form_are_interchangeable(real_crosswalk):
    """Regression: the RNA-seq harmonizer labels a `gene_id` column as
    `ncbi_gene_id`, but published tables put the LOC form there. Searching only
    the numeric column left 90% of a real study unresolved."""
    for id_type in ("ncbi_gene_id", "locus_id"):
        assert resolve_identifier("LOC105320749", id_type).feature_id_standardized == "105320749"
        assert resolve_identifier("105320749", id_type).feature_id_standardized == "105320749"


# --------------------------------------------------------------------------- #
# Crosswalk selection is driven by the study, not by a global setting
# --------------------------------------------------------------------------- #


def test_simulated_study_always_uses_the_demo_crosswalk(monkeypatch):
    monkeypatch.setenv("AREE_CROSSWALK", str(REAL_CROSSWALK))
    set_crosswalk_path(None)
    try:
        assert crosswalk_for_study(simulated=True) == DEMO_CROSSWALK_PATH
    finally:
        set_crosswalk_path(None)


def test_real_study_refuses_to_fall_back_to_the_demo_crosswalk(monkeypatch):
    monkeypatch.delenv("AREE_CROSSWALK", raising=False)
    set_crosswalk_path(None)
    with pytest.raises(FileNotFoundError, match="not simulated"):
        crosswalk_for_study(simulated=False)


def test_real_study_uses_the_selected_crosswalk(monkeypatch):
    monkeypatch.setenv("AREE_CROSSWALK", str(REAL_CROSSWALK))
    set_crosswalk_path(None)
    try:
        assert crosswalk_for_study(simulated=False) == REAL_CROSSWALK
    finally:
        set_crosswalk_path(None)


# --------------------------------------------------------------------------- #
# Real/simulated evidence must never mix
# --------------------------------------------------------------------------- #


def test_evidence_schema_carries_the_simulated_flag():
    from harmonize.schema import EVIDENCE_COLUMNS

    assert "simulated" in EVIDENCE_COLUMNS
    schema = json.loads((REPO_ROOT / "schemas" / "evidence.schema.json").read_text())
    assert schema["properties"]["simulated"]["type"] == "boolean"
    assert "simulated" in schema["required"]


def test_meta_analysis_cannot_pool_across_the_simulated_boundary():
    from meta_analysis.run import GROUP_KEYS

    assert "simulated" in GROUP_KEYS, (
        "without `simulated` in the grouping key, a pooled effect could combine "
        "fabricated demo evidence with evidence from a real study"
    )


# --------------------------------------------------------------------------- #
# Source-table defects
# --------------------------------------------------------------------------- #


def test_converter_drops_a_repeated_header_row():
    """The published sheet '18hrPB+Vc v LOnly' contains a repeated header part-way
    down. It must be dropped and counted, never passed through as an effect size."""
    from intake.supplementary_table import convert_de_table

    frame = pd.DataFrame(
        {
            "Gene ID": ["gene-LOC1|LOC1", "Gene ID", "gene-LOC2|LOC2", None],
            "log2FoldChange": [2.5, "log2FoldChange", -3.0, None],
            "padj": [0.01, "padj", 0.02, None],
        }
    )
    out = Path(__file__).parent / "_tmp_conv.tsv"
    try:
        rep = convert_de_table(
            frame,
            {"gene_id": "Gene ID", "log2FoldChange": "log2FoldChange", "padj": "padj"},
            out_path=out,
        )
        assert rep["rows_written"] == 2
        assert rep["rows_dropped_non_numeric"] == 1
        assert rep["rows_dropped_missing_id_or_effect"] == 1
        assert rep["columns_absent_from_source"] == ["lfcSE", "pvalue"]
        written = pd.read_csv(out, sep="\t")
        assert pd.api.types.is_numeric_dtype(written["log2FoldChange"])
    finally:
        out.unlink(missing_ok=True)


def test_converter_does_not_invent_missing_statistics():
    """A source reporting only padj must yield empty lfcSE/pvalue, not zeros."""
    from intake.supplementary_table import convert_de_table

    frame = pd.DataFrame({"Gene ID": ["gene-LOC1|LOC1"], "log2FoldChange": [2.0], "padj": [0.01]})
    out = Path(__file__).parent / "_tmp_conv2.tsv"
    try:
        convert_de_table(
            frame,
            {"gene_id": "Gene ID", "log2FoldChange": "log2FoldChange", "padj": "padj"},
            out_path=out,
        )
        written = pd.read_csv(out, sep="\t")
        assert written["lfcSE"].isna().all()
        assert written["pvalue"].isna().all()
    finally:
        out.unlink(missing_ok=True)
