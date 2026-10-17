"""Orthogroup-level (cross-species) meta-analysis, WP6 (src/meta_analysis/orthogroups.py)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from conftest import harmonize_all_demo_studies

from meta_analysis.orthogroups import run_orthogroup_meta_analysis
from meta_analysis.pooling import dersimonian_laird

OYSTER, MUSSEL, CLAM = 29159, 29158, 129788


def _gene(feature, taxid, effect, se, k=1, phenotype="thermal_tolerance", simulated=False, p=0.01):
    return {"feature_id_standardized": feature, "phenotype": phenotype, "feature_type": "gene",
            "simulated": simulated, "species_taxid": taxid, "k_studies": k,
            "studies": "|".join(f"S{feature}{i}" for i in range(k)),
            "pooled_effect": effect, "pooled_se": se, "p_value": p, "adjusted_p_value": p}


def _evidence(pairs, simulated=False):
    """Minimal evidence rows carrying (feature, taxid, orthogroup)."""
    return pd.DataFrame([{"feature_id_standardized": f, "simulated": str(simulated), "species_taxid": str(t),
                          "orthogroup_id": og} for f, t, og in pairs])


CLASSES = {"OG1": "single_copy_all", "OG2": "multi_copy", "OG3": "single_copy"}


def test_species_estimates_pool_by_random_effects_across_species():
    genes = pd.DataFrame([_gene("101", OYSTER, 1.0, 0.2, k=2), _gene("201", MUSSEL, 0.6, 0.3, k=2),
                          _gene("301", CLAM, 1.4, 0.25, k=1)])
    ev = _evidence([("101", OYSTER, "OG1"), ("201", MUSSEL, "OG1"), ("301", CLAM, "OG1")])
    result, species, summary = run_orthogroup_meta_analysis(genes, ev, CLASSES)

    row = result.iloc[0]
    expected = dersimonian_laird(np.array([0.6, 1.0, 1.4]), np.array([0.3, 0.2, 0.25]))
    assert row["n_species"] == 3 and row["copy_class"] == "single_copy_all"
    assert row["pooled_effect"] == pytest.approx(expected.pooled_effect)
    assert row["i_squared_between_species"] == pytest.approx(expected.i_squared)
    assert row["direction_consistency_species"] == 1.0
    assert row["k_studies_by_species"] == f"{MUSSEL}:2|{OYSTER}:2|{CLAM}:1"
    assert row["species_with_k2"] == 2 and bool(row["meets_replication_gate"]) is True
    assert len(species) == 3 and summary["pooled"] == 1


def test_one_study_per_species_does_not_meet_the_replication_gate():
    genes = pd.DataFrame([_gene("101", OYSTER, 1.0, 0.2, k=2), _gene("201", MUSSEL, 0.6, 0.3, k=1)])
    ev = _evidence([("101", OYSTER, "OG1"), ("201", MUSSEL, "OG1")])
    result, _, _ = run_orthogroup_meta_analysis(genes, ev, CLASSES)
    assert bool(result.iloc[0]["meets_replication_gate"]) is False


def test_opposite_signs_across_species_are_reported_not_hidden():
    genes = pd.DataFrame([_gene("101", OYSTER, 1.0, 0.2, k=2), _gene("201", MUSSEL, -1.0, 0.2, k=2)])
    ev = _evidence([("101", OYSTER, "OG1"), ("201", MUSSEL, "OG1")])
    row = run_orthogroup_meta_analysis(genes, ev, CLASSES)[0].iloc[0]
    assert row["direction_consistency_species"] == 0.5
    assert row["i_squared_between_species"] > 50


def test_one_species_alone_is_not_a_cross_species_result():
    genes = pd.DataFrame([_gene("101", OYSTER, 1.0, 0.2, k=3)])
    result, _, summary = run_orthogroup_meta_analysis(genes, _evidence([("101", OYSTER, "OG1")]), CLASSES)
    assert len(result) == 0 and summary["single_species"] == 1


def test_multi_copy_orthogroups_are_set_aside_by_default():
    genes = pd.DataFrame([_gene("102", OYSTER, 1.0, 0.2), _gene("202", MUSSEL, 0.8, 0.2)])
    ev = _evidence([("102", OYSTER, "OG2"), ("202", MUSSEL, "OG2")])
    result, _, summary = run_orthogroup_meta_analysis(genes, ev, CLASSES)
    assert len(result) == 0 and summary["multi_copy_excluded"] == 1
    # Opt in: allowed, because each species has evidence for exactly one gene, and labelled.
    result, _, _ = run_orthogroup_meta_analysis(genes, ev, CLASSES, include_multi_copy=True)
    assert result.iloc[0]["copy_class"] == "multi_copy"


def test_paralogs_in_one_species_never_pool():
    genes = pd.DataFrame([_gene("102", OYSTER, 1.0, 0.2), _gene("103", OYSTER, -0.5, 0.2),
                          _gene("202", MUSSEL, 0.8, 0.2)])
    ev = _evidence([("102", OYSTER, "OG2"), ("103", OYSTER, "OG2"), ("202", MUSSEL, "OG2")])
    result, _, summary = run_orthogroup_meta_analysis(genes, ev, CLASSES, include_multi_copy=True)
    assert len(result) == 0 and summary["paralogs_in_species_excluded"] == 1


def test_real_orthogroups_absent_from_the_active_run_are_set_aside():
    genes = pd.DataFrame([_gene("104", OYSTER, 1.0, 0.2), _gene("204", MUSSEL, 0.8, 0.2)])
    ev = _evidence([("104", OYSTER, "OG_OLD"), ("204", MUSSEL, "OG_OLD")])
    result, _, summary = run_orthogroup_meta_analysis(genes, ev, CLASSES)
    assert len(result) == 0 and summary["not_in_active_run_excluded"] == 1


def test_simulated_orthogroups_never_take_a_real_copy_class():
    """Regression: demo IDs (OG0000001...) collide with real OrthoFinder IDs, and
    once took the real run's multi_copy class."""
    genes = pd.DataFrame([_gene("LOCA", OYSTER, 1.0, 0.2, simulated=True),
                          _gene("LOCA", MUSSEL, 0.8, 0.2, simulated=True)])
    ev = _evidence([("LOCA", OYSTER, "OG2"), ("LOCA", MUSSEL, "OG2")], simulated=True)
    result, _, _ = run_orthogroup_meta_analysis(genes, ev, CLASSES)
    assert result.iloc[0]["copy_class"] == "demo_crosswalk"


def test_a_gene_with_two_orthogroups_in_the_evidence_is_an_error():
    genes = pd.DataFrame([_gene("101", OYSTER, 1.0, 0.2), _gene("201", MUSSEL, 0.6, 0.3)])
    ev = _evidence([("101", OYSTER, "OG1"), ("101", OYSTER, "OG3"), ("201", MUSSEL, "OG1")])
    with pytest.raises(ValueError, match="two orthogroup IDs"):
        run_orthogroup_meta_analysis(genes, ev, CLASSES)


def test_adjusted_p_is_within_phenotype_families():
    genes = pd.DataFrame([
        _gene("101", OYSTER, 1.0, 0.2, p=0.01), _gene("201", MUSSEL, 1.0, 0.2, p=0.01),
        _gene("103", OYSTER, 0.2, 0.2, phenotype="survival"), _gene("203", MUSSEL, 0.1, 0.2, phenotype="survival"),
    ])
    ev = _evidence([("101", OYSTER, "OG1"), ("201", MUSSEL, "OG1"), ("103", OYSTER, "OG3"), ("203", MUSSEL, "OG3")])
    result, _, _ = run_orthogroup_meta_analysis(genes, ev, CLASSES)
    assert set(result["n_tests_in_family"]) == {1}


def test_demo_pools_oyster_and_mussel_through_demo_orthogroups(isolated_reports):
    from click.testing import CliRunner

    from aree.cli import main

    harmonize_all_demo_studies()
    out = CliRunner().invoke(main, ["meta-analyze", "--level", "orthogroup", "--feature-type", "gene"])
    assert out.exit_code == 0, out.output
    result = pd.read_csv(isolated_reports["meta_dir"] / "all_phenotypes_gene_orthogroup_meta_analysis.tsv",
                         sep="\t", dtype={"species_taxids": str})
    species = pd.read_csv(isolated_reports["meta_dir"] / "all_phenotypes_gene_orthogroup_species_estimates.tsv",
                          sep="\t")
    assert len(result) > 0
    assert set(result["copy_class"]) == {"demo_crosswalk"}
    assert set(result["species_taxids"]) == {f"{MUSSEL}|{OYSTER}"}
    hsp70 = result[result["orthogroup_id"] == "OG0000001"].iloc[0]
    assert hsp70["phenotype"] == "larval_viability"
    assert hsp70["k_studies_by_species"] == f"{MUSSEL}:1|{OYSTER}:2"
    # The mussel contributes one study, so the demo never clears the gate.
    assert not result["meets_replication_gate"].any()
    # Each species' estimate is exactly its gene-level pooled estimate.
    mussel = species[(species["orthogroup_id"] == "OG0000001") & (species["species_taxid"] == MUSSEL)].iloc[0]
    assert mussel["pooled_effect"] == pytest.approx(2.10)
