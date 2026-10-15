"""Second-species demo (WP4 of docs/mussel_clam_integration_plan.md).

MYTGAL_HEAT01 is a simulated Mytilus galloprovincialis larval study that
reuses the oyster demo's gene identifiers under the same phenotype
(larval_viability) as the oyster studies GIGAS_OA02 and GIGAS_LARV05. These
tests pin down what that must and must not change:

* it forms its own species partition in evidence, meta-analysis, candidates and
  evidence cards;
* every oyster result, down to adjusted p-values and candidate scores, is
  identical with and without it. A score or multiple-testing family that
  silently spanned species would fail here.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from conftest import ALL_DEMO_STUDY_IDS, OYSTER_TAXID, oyster

from harmonize.core import harmonize_study, load_evidence_table
from meta_analysis.run import run_meta_analysis
from prioritize.rank import build_candidates
from reporting.evidence_cards import build_evidence_cards

MUSSEL_STUDY = "MYTGAL_HEAT01"
MUSSEL_TAXID = 29158
HSP70, SOD1 = "LOC105333935", "LOC105331241"


def _harmonize(study_ids):
    for sid in study_ids:
        harmonize_study(sid, date_generated="2026-01-01")


def _mussel(df):
    return df[df["species_taxid"].astype(int) == MUSSEL_TAXID]


def _sorted(df):
    keys = ["feature_id_standardized", "phenotype", "feature_type"]
    return df.sort_values(keys).reset_index(drop=True)


def test_mussel_demo_study_is_a_registered_simulated_study():
    from intake.schema_validate import validate_study_file

    assert MUSSEL_STUDY in ALL_DEMO_STUDY_IDS
    result = validate_study_file(Path("registry/studies") / f"{MUSSEL_STUDY}.yaml")
    assert result.valid, result.errors


def test_mussel_evidence_is_its_own_species_and_carries_the_unverified_flag(isolated_reports):
    _harmonize(ALL_DEMO_STUDY_IDS)
    evidence = load_evidence_table()
    mussel = evidence[evidence["study_id"] == MUSSEL_STUDY]

    assert len(mussel) == 10
    assert set(mussel["species_taxid"].astype(int)) == {MUSSEL_TAXID}
    assert set(mussel["species"]) == {"Mytilus galloprovincialis"}
    assert (mussel["mapping_confidence"] == "exact").all()
    assert mussel["quality_flags"].str.contains("species_assignment_unverified").all()
    # Only the mussel study contributes mussel evidence.
    assert set(_mussel(evidence)["study_id"]) == {MUSSEL_STUDY}


def test_shared_identifiers_form_separate_pools_per_species(isolated_reports):
    _harmonize(ALL_DEMO_STUDY_IDS)
    meta = run_meta_analysis(phenotype="larval_viability", feature_type="gene")
    hsp70 = meta[meta["feature_id_standardized"] == HSP70]

    assert sorted(hsp70["species_taxid"].astype(int)) == [MUSSEL_TAXID, OYSTER_TAXID]
    mussel = _mussel(hsp70).iloc[0]
    assert mussel["k_studies"] == 1
    assert mussel["pooled_effect"] == pytest.approx(2.10)
    assert oyster(hsp70).iloc[0]["k_studies"] == 2


def test_oyster_meta_analysis_is_identical_with_and_without_the_mussel_study(isolated_reports):
    oyster_only = [s for s in ALL_DEMO_STUDY_IDS if s != MUSSEL_STUDY]

    _harmonize(oyster_only)
    before = _sorted(oyster(run_meta_analysis(feature_type="gene")))
    _harmonize([MUSSEL_STUDY])
    after_all = run_meta_analysis(feature_type="gene")
    after = _sorted(oyster(after_all))

    pd.testing.assert_frame_equal(before, after)
    # The mussel study added rows, it did not merge into existing ones.
    assert len(after_all) == len(before) + len(_mussel(after_all))


def test_oyster_candidates_and_scores_are_identical_with_and_without_the_mussel_study(isolated_reports):
    oyster_only = [s for s in ALL_DEMO_STUDY_IDS if s != MUSSEL_STUDY]

    _harmonize(oyster_only)
    before = _sorted(oyster(build_candidates(
        run_meta_analysis(feature_type="gene"), load_evidence_table())))
    _harmonize([MUSSEL_STUDY])
    after = _sorted(oyster(build_candidates(
        run_meta_analysis(feature_type="gene"), load_evidence_table())))

    pd.testing.assert_frame_equal(before, after)


def test_single_study_mussel_candidates_are_never_high_priority(isolated_reports):
    _harmonize(ALL_DEMO_STUDY_IDS)
    candidates = build_candidates(run_meta_analysis(feature_type="gene"), load_evidence_table())
    mussel = _mussel(candidates)

    assert len(mussel) > 0
    assert (mussel["k_studies"] == 1).all()
    assert not (mussel["tier"] == "high_priority_cross_study").any()


def test_mussel_and_oyster_get_separate_cards_that_name_their_species(isolated_reports):
    _harmonize(ALL_DEMO_STUDY_IDS)
    index = build_evidence_cards(phenotype="larval_viability", feature_type="gene")
    sod1 = {int(e["species_taxid"]): e for e in index if e["feature_id_standardized"] == SOD1}

    assert set(sod1) == {OYSTER_TAXID, MUSSEL_TAXID}
    assert sod1[OYSTER_TAXID]["card_file"] != sod1[MUSSEL_TAXID]["card_file"]

    mussel_card = Path(sod1[MUSSEL_TAXID]["card_file"]).read_text()
    oyster_card = Path(sod1[OYSTER_TAXID]["card_file"]).read_text()
    assert "*Mytilus galloprovincialis* (Mediterranean mussel)" in mussel_card
    assert "*Magallana gigas* (Pacific oyster)" in oyster_card
    assert MUSSEL_STUDY in mussel_card and MUSSEL_STUDY not in oyster_card
    assert "GIGAS_OA02" not in mussel_card
    # The oyster studies disagree on sod1; the mussel's opposite sign must not
    # appear as a third, conflicting study on the oyster card.
    assert "Conflicting direction" in oyster_card
    for card in (mussel_card, oyster_card):
        assert "never pooled with this one" in card
