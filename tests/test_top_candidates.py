from conftest import harmonize_all_demo_studies

from reporting.evidence_cards import build_evidence_cards_report
from reporting.top_candidates import build_top_candidates_summary, top_candidates_by_phenotype


def test_top_candidates_by_phenotype_respects_n_and_ranking(isolated_reports):
    harmonize_all_demo_studies()
    report = build_evidence_cards_report(phenotype="larval_viability", feature_type="gene", all_cards=True)
    import pandas as pd
    candidates = pd.read_csv(report.candidates_path, sep="\t")

    by_phenotype = top_candidates_by_phenotype(candidates, n=2)
    group = by_phenotype["larval_viability"]
    assert len(group) <= 2
    assert list(group["score"]) == sorted(group["score"], reverse=True)


def test_build_top_candidates_summary_writes_markdown(isolated_reports, tmp_path):
    harmonize_all_demo_studies()
    report = build_evidence_cards_report(phenotype="larval_viability", feature_type="gene", all_cards=True)
    out_path = tmp_path / "summary.md"

    result_path = build_top_candidates_summary(n=3, candidates_path=report.candidates_path, out_path=out_path)

    assert result_path == out_path
    content = out_path.read_text()
    assert "Top candidate biomarkers by phenotype" in content
    assert "larval_viability" in content
    assert "not validation" in content
