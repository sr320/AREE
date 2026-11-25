"""Unit tests for harmonization dispatch helpers."""
from __future__ import annotations

import pytest

from common import DEMO_STUDIES_DIR, REPO_ROOT, load_yaml
from harmonize.core import _comparison_assay_type, harmonize_batch, load_evidence_table


def test_comparison_assay_type_single_and_multi_assay():
    single = {"study_id": "S1", "assay_type": ["rnaseq"]}
    assert _comparison_assay_type(single, {"comparison_id": "c1"}) == "rnaseq"

    multi = {"study_id": "S2", "assay_type": ["rnaseq", "proteomics"]}
    assert _comparison_assay_type(multi, {"comparison_id": "c1", "assay_type": "proteomics"}) == "proteomics"

    with pytest.raises(ValueError, match="multiple assay_types"):
        _comparison_assay_type(multi, {"comparison_id": "c1"})

    with pytest.raises(ValueError, match="not in study assay_type"):
        _comparison_assay_type(multi, {"comparison_id": "c1", "assay_type": "metabolomics"})


def test_harmonize_batch_writes_all_comparisons_in_one_upsert(isolated_reports):
    jobs = []
    for study_id in ("GIGAS_HEAT01", "GIGAS_OA02"):
        study = load_yaml(DEMO_STUDIES_DIR / f"{study_id}.yaml")
        for comp in study["comparisons"]:
            jobs.append((study_id, comp["comparison_id"], REPO_ROOT / comp["results_file"]))

    counts = harmonize_batch(jobs, date_generated="2026-01-01")
    assert len(counts) == len(jobs)
    evidence = load_evidence_table()
    assert len(evidence) == sum(n for _, _, n in counts)
    assert set(evidence["study_id"]) == {"GIGAS_HEAT01", "GIGAS_OA02"}


def test_use_crosswalk_restores_previous_state(tmp_path):
    from harmonize.identifiers import (
        DEMO_CROSSWALK_PATH,
        active_crosswalk_path,
        set_crosswalk_path,
        use_crosswalk,
    )

    alt = tmp_path / "alt_crosswalk.tsv"
    alt.write_text("ncbi_gene_id\tgene_symbol\n1\tG1\n")
    set_crosswalk_path(None)
    assert active_crosswalk_path() == DEMO_CROSSWALK_PATH
    with use_crosswalk(alt):
        assert active_crosswalk_path() == alt
    assert active_crosswalk_path() == DEMO_CROSSWALK_PATH


def test_statistical_guards_and_deterministic_sort():
    import pandas as pd

    from meta_analysis.pooling import dersimonian_laird
    from prioritize.rank import sort_candidates
    from prioritize.scoring import _clip01, _saturating

    with pytest.raises(ValueError, match="strictly positive"):
        dersimonian_laird([1.0], [0.0])
    with pytest.raises(ValueError, match="finite"):
        dersimonian_laird([float("nan")], [0.2])
    with pytest.raises(ValueError, match="finite"):
        _clip01(float("nan"))
    with pytest.raises(ValueError, match="finite"):
        _saturating(float("nan"), 2.0)

    df = pd.DataFrame({
        "tier": ["emerging", "emerging", "emerging"],
        "score": [50.0, 50.0, 50.0],
        "adjusted_p_value": [0.04, 0.01, 0.01],
        "feature_id_standardized": ["G3", "G2", "G1"],
    })
    assert sort_candidates(df)["feature_id_standardized"].tolist() == ["G1", "G2", "G3"]
