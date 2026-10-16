"""The methylation harmonizer carries a DMR table's standard error, p-value and
method when the DMR step supplies them, and leaves older tables unchanged.

Without a standard error and an unadjusted p-value a methylation result
cannot enter inverse-variance pooling. methylKit tables carry neither;
DMR_POOLED_TILES (dmr_engine limma or dss) writes both.
"""
from __future__ import annotations

import pandas as pd

from common import REPO_ROOT, STUDIES_DIR, load_yaml

STUDY_ID = "GIGAS_PATH03"


def _demo_table():
    study = load_yaml(STUDIES_DIR / f"{STUDY_ID}.yaml")
    (comp,) = study["comparisons"]
    return comp["comparison_id"], pd.read_csv(REPO_ROOT / comp["results_file"], sep="\t")


def test_methylkit_table_stays_unpoolable_and_labelled_methylkit(isolated_reports, tmp_path):
    from harmonize.core import harmonize_processed_table

    cid, table = _demo_table()
    path = tmp_path / "plain_dmr_standardized.tsv"
    table.to_csv(path, sep="\t", index=False)
    ev = harmonize_processed_table(STUDY_ID, path, date_generated="2026-01-01", comparison_id=cid)
    assert ev["standard_error"].isna().all()
    assert ev["p_value"].isna().all()
    assert set(ev["analysis_method"]) == {"methylKit"}


def test_pooled_tile_columns_reach_the_evidence_table(isolated_reports, tmp_path):
    from harmonize.core import harmonize_processed_table

    cid, table = _demo_table()
    table["pvalue"] = 0.01
    table["se_meth_diff_percent"] = 2.5
    table["analysis_method"] = "limma_pooled_tiles"
    path = tmp_path / "pooled_dmr_standardized.tsv"
    table.to_csv(path, sep="\t", index=False)
    ev = harmonize_processed_table(STUDY_ID, path, date_generated="2026-01-01", comparison_id=cid)
    assert (ev["standard_error"] == 2.5).all()
    assert (ev["p_value"] == 0.01).all()
    assert set(ev["analysis_method"]) == {"limma_pooled_tiles"}
    assert ev["evidence_id"].str.contains("limma_pooled_tiles").all()
