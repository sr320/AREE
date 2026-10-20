"""Curation tests for WANG2026_ERK_HEAT.

Real proteomics study (PRIDE PXD065668 + PXD065706; Wang et al. 2026, Commun
Biol), registered but deliberately not reanalyzed. These tests pin the facts
behind that decision: 4 groups x 3 pools, quantification only in proprietary
Spectronaut archives, supplementary tables that hold significant features
only, and proteome runs acquired in condition blocks.
"""
from __future__ import annotations

import hashlib
import json

import pandas as pd
import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "WANG2026_ERK_HEAT.yaml"
STUDY_DIR = DATA_DIR / "studies" / "WANG2026_ERK_HEAT"


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def provenance():
    return json.loads((STUDY_DIR / "pride_provenance.json").read_text())


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_both_projects_are_four_groups_of_three_pools_with_a_spectronaut_archive(provenance):
    for dataset in ("proteome", "phosphoproteome"):
        p = provenance["pride"][dataset]
        assert p["pools_per_group"] == {"ERK_0": 3, "ERK_12": 3, "SW_0": 3, "SW_12": 3}
        assert p["files_by_category"] == {"RAW": 12, "SEARCH": 1}
        assert len(p["search_files"]) == 1 and p["search_files"][0].endswith(".sne")


def test_supplementary_tables_hold_significant_features_only(provenance):
    for dataset in ("proteome", "phosphoproteome"):
        for sheet in provenance["supplementary_tables"][dataset]["sheets"].values():
            assert {r.lower() for r in sheet["regulation"]} <= {"up", "down"}
    sheets = provenance["supplementary_tables"]["proteome"]["sheets"]
    assert sheets["Sheet 1"]["n_rows"] == 476 and sheets["Sheet 3"]["n_rows"] == 348
    assert "CC BY-NC-ND" in provenance["paper_license"]
    assert not list(STUDY_DIR.glob("*.xlsx"))


def test_proteome_runs_were_acquired_in_condition_blocks(provenance):
    order = [r.rsplit("_", 1)[0] for r in provenance["proteome_acquisition_order"]]
    assert order == ["seawater_0h"] * 3 + ["seawater_12h"] * 3 + ["SCH772984_0h"] * 3 + ["SCH772984_12h"] * 3
    runs = pd.read_csv(STUDY_DIR / "run_order.tsv", sep="\t")
    assert len(runs) == 12


def test_registered_not_reanalyzed(study):
    assert study["qc_status"] == "not_started" and study["analysis_status"] == "not_started"
    assert all(c["results_file"] is None for c in study["comparisons"])
    assert "batch_confounded" in study["quality_flags"]
    assert study["accessions"]["proteomexchange"] == "PXD065668; PXD065706"


def test_comparisons(study):
    comps = {c["comparison_id"]: c for c in study["comparisons"]}
    assert [cid for cid, c in comps.items() if c.get("meta_analysis_primary")] == ["heat_12h_vs_0h_seawater"]
    assert comps["heat_12h_vs_0h_seawater"]["resilience_classification"] == "exposure_only"
    inhib = comps["erk_inhibitor_vs_seawater_at_12h_heat"]
    assert (inhib["phenotype"], inhib["phenotype_direction"]) == ("thermal_tolerance", "decrease")
    for c in comps.values():
        assert (c["sample_size"], c["biological_replicates"]) == (6, 3)


def test_run_order_checksum_matches_provenance(provenance):
    for rel, digest in provenance["outputs_sha256"].items():
        assert hashlib.sha256((REPO_ROOT / rel).read_bytes()).hexdigest() == digest
