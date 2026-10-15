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



def test_summary_names_the_species_of_every_row(isolated_reports, tmp_path):
    """The demo holds two species that share feature IDs, so the same ID can rank
    twice under one phenotype; each row must say which species it is."""
    harmonize_all_demo_studies()
    report = build_evidence_cards_report(phenotype="larval_viability", feature_type="gene", all_cards=True)
    out = build_top_candidates_summary(n=50, candidates_path=report.candidates_path,
                                       out_path=tmp_path / "summary.md")

    hsp70_rows = [line for line in out.read_text().splitlines() if "| LOC105333935 " in line]
    assert len(hsp70_rows) == 2
    assert sum("Mytilus galloprovincialis" in r for r in hsp70_rows) == 1
    assert sum("Magallana gigas" in r for r in hsp70_rows) == 1
