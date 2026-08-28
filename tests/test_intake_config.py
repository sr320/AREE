"""Tests for the declarative intake step (published artifact -> AREE result files).

The headline guarantee is `test_committed_real_study_reproduces_from_source`:
the derived TSVs committed for the real study must still follow, byte for byte,
from the committed published artifact. Everything else here checks that a
malformed or drifted intake fails loudly instead of quietly degrading.
"""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import pytest
import yaml

from common import REPO_ROOT
from intake.run_intake import IntakeError, run_intake

HESSER_CONFIG = REPO_ROOT / "data/studies/HESSER2024_VCOR/intake.yaml"
HESSER_DIR = REPO_ROOT / "data/studies/HESSER2024_VCOR"


def _tree_fingerprint(directory: Path) -> dict:
    """sha256 of every file under `directory`, for asserting nothing was touched."""
    return {
        str(p.relative_to(directory)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(directory.rglob("*"))
        if p.is_file()
    }


@pytest.fixture
def hesser_cfg():
    pytest.importorskip("xlrd", reason="spreadsheet intake needs: pip install -e '.[intake]'")
    return yaml.safe_load(HESSER_CONFIG.read_text())


def _write_cfg(tmp_path: Path, cfg: dict) -> Path:
    path = tmp_path / "intake.yaml"
    path.write_text(yaml.safe_dump(cfg))
    return path


# --------------------------------------------------------------------------- #
# The reproducibility guarantee
# --------------------------------------------------------------------------- #


def test_committed_real_study_reproduces_from_source(hesser_cfg):
    """The committed HESSER2024_VCOR result files must regenerate exactly."""
    report = run_intake(HESSER_CONFIG, check=True)
    assert report["mismatches"] == [], (
        "the committed derived files no longer follow from the committed source artifact"
    )
    assert len(report["conversions"]) == 2


def test_check_mode_does_not_modify_the_repository(hesser_cfg):
    before = _tree_fingerprint(HESSER_DIR)
    run_intake(HESSER_CONFIG, check=True)
    assert _tree_fingerprint(HESSER_DIR) == before


def test_check_reports_row_counts_matching_committed_provenance(hesser_cfg):
    by_name = {
        Path(c["output_file"]).name: c for c in run_intake(HESSER_CONFIG, check=True)["conversions"]
    }
    assert by_name["HESSER2024_VCOR_vcor_vs_larvae_only_dge.tsv"]["rows_written"] == 129
    probiotic = by_name["HESSER2024_VCOR_probiotic_plus_vcor_vs_larvae_only_dge.tsv"]
    assert probiotic["rows_written"] == 222
    # The repeated header row part-way down the published sheet.
    assert probiotic["rows_dropped_non_numeric"] == 1
    # The source reports neither, and intake must never invent them.
    assert probiotic["columns_absent_from_source"] == ["lfcSE", "pvalue"]


# --------------------------------------------------------------------------- #
# Drift and malformed configs must fail loudly
# --------------------------------------------------------------------------- #


def test_source_checksum_mismatch_is_fatal(hesser_cfg, tmp_path):
    hesser_cfg["source"]["sha256"] = "0" * 64
    with pytest.raises(IntakeError, match="checksum mismatch"):
        run_intake(_write_cfg(tmp_path, hesser_cfg), check=True)


def test_missing_source_artifact_is_reported_with_its_url(hesser_cfg, tmp_path):
    hesser_cfg["source"]["local_copy"] = "data/studies/NOPE/missing.xls"
    with pytest.raises(IntakeError, match="Source artifact not found"):
        run_intake(_write_cfg(tmp_path, hesser_cfg), check=True)


def test_transformation_note_parsed_as_a_mapping_is_rejected(hesser_cfg, tmp_path):
    """An unquoted note containing ': ' becomes a YAML dict; that must not reach provenance."""
    hesser_cfg["transformation_notes"] = [{"Sheet X was not converted": "it has no stressor"}]
    with pytest.raises(IntakeError, match="must be strings"):
        run_intake(_write_cfg(tmp_path, hesser_cfg), check=True)


def test_unknown_aree_column_in_column_map_is_rejected(hesser_cfg, tmp_path):
    hesser_cfg["conversions"][0]["column_map"]["log2FC"] = "log2FoldChange"
    with pytest.raises(IntakeError, match="unknown AREE column"):
        run_intake(_write_cfg(tmp_path, hesser_cfg), check=True)


def test_conversion_without_a_gene_id_mapping_is_rejected(hesser_cfg, tmp_path):
    del hesser_cfg["conversions"][0]["column_map"]["gene_id"]
    with pytest.raises(IntakeError, match="must map 'gene_id'"):
        run_intake(_write_cfg(tmp_path, hesser_cfg), check=True)


def test_missing_required_top_level_key_is_rejected(hesser_cfg, tmp_path):
    del hesser_cfg["conversions"]
    with pytest.raises(IntakeError, match="missing required key 'conversions'"):
        run_intake(_write_cfg(tmp_path, hesser_cfg), check=True)


def test_nonexistent_sheet_names_the_sheet(hesser_cfg, tmp_path):
    hesser_cfg["conversions"][0]["source_sheet"] = "No Such Sheet"
    with pytest.raises(IntakeError, match="No Such Sheet"):
        run_intake(_write_cfg(tmp_path, hesser_cfg), check=True)


# --------------------------------------------------------------------------- #
# Non-spreadsheet sources
# --------------------------------------------------------------------------- #


def test_csv_source_needs_no_spreadsheet_engine(tmp_path):
    """A published CSV must work without the optional [intake] extra."""
    source = tmp_path / "supp.csv"
    source.write_text("Gene ID,log2FoldChange,padj\ngene-LOC1|LOC1,2.5,0.01\n")
    cfg = {
        "study_id": "CSV_DEMO",
        "source": {"local_copy": str(source)},
        "output_dir": str(tmp_path / "out"),
        "provenance_file": str(tmp_path / "out" / "intake_provenance.json"),
        "conversions": [
            {
                "output_file": "CSV_DEMO_dge.tsv",
                "column_map": {"gene_id": "Gene ID", "log2FoldChange": "log2FoldChange", "padj": "padj"},
            }
        ],
    }
    report = run_intake(_write_cfg(tmp_path, cfg))
    assert report["conversions"][0]["rows_written"] == 1
    assert (tmp_path / "out" / "CSV_DEMO_dge.tsv").exists()


def test_hand_edited_result_file_is_detected(hesser_cfg, tmp_path):
    """If someone edits a derived TSV by hand, --check must catch it."""
    work = tmp_path / "HESSER2024_VCOR"
    shutil.copytree(HESSER_DIR, work)
    edited = work / "HESSER2024_VCOR_vcor_vs_larvae_only_dge.tsv"
    edited.write_text(edited.read_text() + "gene-LOC999|LOC999\t9.9\t\t\t0.001\n")

    hesser_cfg["source"]["local_copy"] = str(work / "_source" / "fimmu-15-1380089-Table_S4.xls")
    hesser_cfg["output_dir"] = str(work)
    hesser_cfg["provenance_file"] = str(work / "intake_provenance.json")

    report = run_intake(_write_cfg(tmp_path, hesser_cfg), check=True)
    assert any("does not match what the config regenerates" in m for m in report["mismatches"])


def test_stale_provenance_checksum_is_detected(hesser_cfg, tmp_path):
    """Provenance that no longer describes the committed source must not pass silently."""
    work = tmp_path / "HESSER2024_VCOR"
    shutil.copytree(HESSER_DIR, work)
    prov_path = work / "intake_provenance.json"
    prov = json.loads(prov_path.read_text())
    prov["conversions"][0]["output_sha256"] = "0" * 64
    prov_path.write_text(json.dumps(prov, indent=2))

    hesser_cfg["source"]["local_copy"] = str(work / "_source" / "fimmu-15-1380089-Table_S4.xls")
    hesser_cfg["output_dir"] = str(work)
    hesser_cfg["provenance_file"] = str(prov_path)

    report = run_intake(_write_cfg(tmp_path, hesser_cfg), check=True)
    assert any("!= provenance" in m for m in report["mismatches"])
