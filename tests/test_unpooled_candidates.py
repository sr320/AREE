"""Evidence that cannot be pooled (no standard error, no unadjusted p) must stay
out of every pooled estimate yet remain visible as an unscored candidate.

HESSER2024_VCOR publishes only a significance-filtered DEG list with log2FC and
padj; before this, its `disease_resistance` comparison — the only real
resilience-classified comparison in the registry — vanished from the ranking."""
import math
from pathlib import Path

import pandas as pd
from conftest import harmonize_all_demo_studies

from harmonize.core import load_evidence_table
from meta_analysis.run import run_meta_analysis
from reporting.evidence_cards import build_evidence_cards_report


def _write_evidence(isolated_reports, rows):
    path = isolated_reports["evidence_table_path"]
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def _rows(source, specs):
    rows = []
    for evidence_id, study_id, phenotype, se, p, padj in specs:
        row = source.copy()
        row["evidence_id"] = evidence_id
        row["study_id"] = study_id
        row["comparison_id"] = f"{study_id}_comparison"
        row["phenotype"] = phenotype
        row["standard_error"] = se
        row["p_value"] = p
        row["adjusted_p_value"] = padj
        rows.append(row)
    return rows


def _source():
    evidence = load_evidence_table()
    return evidence[
        (evidence["feature_id_standardized"] == "LOC105331241")
        & (evidence["phenotype"] == "larval_viability")
    ].iloc[0]


def test_padj_only_record_is_never_pooled(isolated_reports):
    harmonize_all_demo_studies()
    _write_evidence(isolated_reports, _rows(_source(), [
        ("pooled", "POOLED_STUDY", "larval_viability", 0.2, 0.01, 0.02),
        ("padj_only", "FILTERED_STUDY", "larval_viability", None, None, 1e-12),
    ]))
    result = run_meta_analysis(phenotype="larval_viability", feature_type="gene").iloc[0]
    assert result["studies"] == "POOLED_STUDY"
    assert result["n_excluded_unpoolable"] == 1
    assert result["excluded_studies"] == "FILTERED_STUDY"


def test_unpoolable_only_feature_gets_an_unscored_emerging_row_and_card(isolated_reports):
    harmonize_all_demo_studies()
    _write_evidence(isolated_reports, _rows(_source(), [
        ("pooled", "POOLED_STUDY", "larval_viability", 0.2, 0.01, 0.02),
        ("padj_only", "FILTERED_STUDY", "disease_resistance", None, None, 1e-12),
    ]))
    report = build_evidence_cards_report()
    candidates = pd.read_csv(report.candidates_path, sep="\t")

    unpooled = candidates[candidates["phenotype"] == "disease_resistance"]
    assert len(unpooled) == 1
    row = unpooled.iloc[0]
    assert row["pooling_status"] == "unpooled"
    assert row["tier"] == "emerging"
    assert math.isnan(row["score"])
    assert row["k_studies"] == 0
    assert row["excluded_studies"] == "FILTERED_STUDY"

    pooled = candidates[candidates["phenotype"] == "larval_viability"].iloc[0]
    assert pooled["pooling_status"] == "pooled"
    # Unscored rows sort after every scored candidate.
    assert candidates.index[candidates["pooling_status"] == "unpooled"].min() > \
        candidates.index[candidates["pooling_status"] == "pooled"].max()

    card = Path(row["card_file"]).read_text()
    assert "not scored (unpooled)" in card
    assert "significance-filtered" in card
    assert "FILTERED_STUDY" in card


def test_feature_with_a_pooled_candidate_gets_no_extra_unpooled_row(isolated_reports):
    harmonize_all_demo_studies()
    _write_evidence(isolated_reports, _rows(_source(), [
        ("pooled", "POOLED_STUDY", "larval_viability", 0.2, 0.01, 0.02),
        ("padj_only", "FILTERED_STUDY", "larval_viability", None, None, 1e-12),
    ]))
    report = build_evidence_cards_report()
    candidates = pd.read_csv(report.candidates_path, sep="\t")
    assert len(candidates) == 1
    assert candidates.iloc[0]["pooling_status"] == "pooled"
    card = Path(candidates.iloc[0]["card_file"]).read_text()
    assert "Records left out of the pool" in card
    assert "FILTERED_STUDY" in card
