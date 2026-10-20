"""contrast_type keeps baseline differences between resilient and susceptible
groups apart from responses to exposure, from the registry through ranking."""
import pandas as pd
import pytest
import yaml
from conftest import harmonize_all_demo_studies

from common import STUDIES_DIR, study_files, study_yaml_path
from harmonize.core import load_evidence_table
from intake.schema_validate import validate_study_file
from meta_analysis.run import require_contrast_type, run_meta_analysis
from prioritize.rank import CONTRAST_TYPES, build_candidates
from reporting.evidence_cards import _card_markdown, _card_slug
from reporting.top_candidates import build_top_candidates_summary

BASELINE = "baseline_phenotype_difference"
RESPONSE = "exposure_response"


def _registered_comparisons():
    for path in study_files():  # real and demo
        study = yaml.safe_load(path.read_text())
        for comp in study["comparisons"]:
            yield study["study_id"], comp


def test_every_registered_comparison_declares_a_contrast_type():
    contrasts = {(sid, c["comparison_id"]): c.get("contrast_type") for sid, c in _registered_comparisons()}
    assert set(contrasts.values()) <= set(CONTRAST_TYPES)
    baseline = {k for k, v in contrasts.items() if v == BASELINE}
    assert ("LI2018_HEAT_TIDAL", "intertidal_vs_subtidal_0h") in baseline
    # Starvation sampled before infection is a treatment, not a baseline difference.
    assert contrasts[("DUPERRET2025_STARVE_RNASEQ", "starved_vs_fed_0hpc")] == RESPONSE


def test_baseline_contrasts_are_classified_as_resilience():
    for sid, comp in _registered_comparisons():
        if comp["contrast_type"] == BASELINE:
            assert comp["resilience_classification"] == "resilience", (sid, comp["comparison_id"])


def _copy_study(tmp_path, study_id, edit):
    study = yaml.safe_load(study_yaml_path(study_id).read_text())
    edit(study)
    path = tmp_path / f"{study_id}.yaml"
    path.write_text(yaml.safe_dump(study, sort_keys=False))
    return path


def test_comparison_without_contrast_type_fails_validation(tmp_path):
    path = _copy_study(tmp_path, "GIGAS_HEAT01", lambda s: s["comparisons"][0].pop("contrast_type"))
    result = validate_study_file(path)
    assert not result.valid
    assert any("contrast_type" in e for e in result.errors)


def test_baseline_contrast_on_a_non_resilience_comparison_warns(tmp_path):
    def edit(study):
        study["comparisons"][0]["contrast_type"] = BASELINE
        study["comparisons"][0]["resilience_classification"] = "stress_response"
    result = validate_study_file(_copy_study(tmp_path, "GIGAS_HEAT01", edit))
    assert any("baseline_phenotype_difference" in w for w in result.warnings)


def test_evidence_without_contrast_type_is_refused():
    with pytest.raises(ValueError, match="re-harmonize"):
        require_contrast_type(pd.DataFrame({"study_id": ["OLD"], "effect_size": [1.0]}))
    with pytest.raises(ValueError, match="OLD"):
        require_contrast_type(pd.DataFrame({"study_id": ["OLD", "NEW"], "contrast_type": [None, RESPONSE]}))


@pytest.fixture
def two_contrast_study(isolated_reports, tmp_path, monkeypatch):
    """One study, one gene, one phenotype: a baseline contrast and two exposure
    responses, one of them flagged primary."""
    import meta_analysis.run as meta_run

    harmonize_all_demo_studies()
    evidence = load_evidence_table()
    source = evidence[
        (evidence["feature_id_standardized"] == "LOC105331241")
        & (evidence["phenotype"] == "larval_viability")
    ].iloc[0]

    rows = []
    for comparison_id, contrast, effect in [("baseline", BASELINE, -0.9), ("response_a", RESPONSE, 1.4),
                                            ("response_b", RESPONSE, 1.1)]:
        row = source.copy()
        row["evidence_id"] = comparison_id
        row["study_id"] = "TWO_CONTRAST"
        row["comparison_id"] = comparison_id
        row["contrast_type"] = contrast
        row["resilience_classification"] = "resilience"
        row["effect_size"] = effect
        row["standard_error"] = 0.2
        row["p_value"] = 1e-4
        row["adjusted_p_value"] = 1e-3
        rows.append(row)
    path = isolated_reports["evidence_table_path"]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)

    studies_dir = tmp_path / "studies"
    studies_dir.mkdir()
    (studies_dir / "TWO_CONTRAST.yaml").write_text(
        "study_id: TWO_CONTRAST\n"
        "comparisons:\n"
        "  - comparison_id: baseline\n"
        f"    contrast_type: {BASELINE}\n"
        "  - comparison_id: response_a\n"
        f"    contrast_type: {RESPONSE}\n"
        "    meta_analysis_primary: true\n"
        "  - comparison_id: response_b\n"
        f"    contrast_type: {RESPONSE}\n"
    )
    monkeypatch.setattr(meta_run, "study_files", lambda: sorted(studies_dir.glob("*.yaml")))
    return load_evidence_table()


def test_contrast_types_pool_apart_and_primary_is_chosen_within_each(two_contrast_study):
    result = run_meta_analysis(phenotype="larval_viability", feature_type="gene").set_index("contrast_type")
    assert set(result.index) == {BASELINE, RESPONSE}
    # The flagged exposure response excludes its sibling, not the baseline contrast.
    assert result.loc[RESPONSE, "contributing_evidence_ids"] == "response_a"
    assert result.loc[RESPONSE, "n_excluded_non_primary"] == 1
    assert result.loc[BASELINE, "contributing_evidence_ids"] == "baseline"
    assert result.loc[BASELINE, "n_excluded_non_primary"] == 0
    # Each contrast type is its own multiple-testing family.
    assert set(result["n_tests_in_family"]) == {1}


def test_baseline_candidate_is_flagged_predictive_and_gets_its_own_card(two_contrast_study, tmp_path):
    meta = run_meta_analysis(phenotype="larval_viability", feature_type="gene")
    candidates = build_candidates(meta, two_contrast_study).set_index("contrast_type")
    baseline, response = candidates.loc[BASELINE].to_dict(), candidates.loc[RESPONSE].to_dict()

    assert baseline["is_predictive_baseline"] and not response["is_predictive_baseline"]
    assert baseline["evidence_class"] == "resilience_associated"
    baseline["contrast_type"], response["contrast_type"] = BASELINE, RESPONSE
    assert _card_slug(baseline).endswith("_baseline_phenotype_difference")
    assert not _card_slug(response).endswith(RESPONSE)  # existing card names unchanged
    assert _card_slug(baseline) != _card_slug(response)
    card = _card_markdown(baseline, [], [], "Larval viability")
    assert "Predictive baseline candidate" in card and "common unstressed conditions" in card

    path = tmp_path / "candidates.tsv"
    candidates.reset_index().assign(card_file="").to_csv(path, sep="\t", index=False)
    # The fixture is built from demo evidence, so it is simulated.
    summary = build_top_candidates_summary(candidates_path=path, out_path=tmp_path / "s.md",
                                           include_simulated=True).read_text()
    section = summary.split("## Predictive baseline candidates", 1)[1].split("\n## ", 1)[0]
    assert "LOC105331241" in section and BASELINE in section and RESPONSE not in section


def test_summary_says_so_when_there_is_no_baseline_candidate(isolated_reports, tmp_path):
    from reporting.evidence_cards import build_evidence_cards_report

    harmonize_all_demo_studies()
    report = build_evidence_cards_report(phenotype="larval_viability", feature_type="gene", all_cards=True)
    summary = build_top_candidates_summary(candidates_path=report.candidates_path,
                                           out_path=tmp_path / "s.md").read_text()
    assert "No candidate rests on a baseline contrast with a resilience outcome yet." in summary


def test_template_documents_the_choice():
    text = (STUDIES_DIR / "_TEMPLATE.yaml").read_text()
    assert "contrast_type:" in text and "docs/contrast_types.md" in text
