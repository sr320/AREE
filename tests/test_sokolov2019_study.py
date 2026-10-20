"""Curation tests for SOKOLOV2019_HYPOX_MITO.

Real proteomics study (PRIDE PXD011365; Sokolov et al. 2019, J Proteomics),
registered but deliberately not reanalyzed. These tests pin the facts behind
that decision: the total proteome is 3 oysters per condition acquired one
oyster at a time in condition blocks, the phosphoproteome has no condition
key, and nothing quantitative beyond Scaffold mzIdentML exports was deposited.
"""
from __future__ import annotations

import hashlib
import json

import pandas as pd
import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "SOKOLOV2019_HYPOX_MITO.yaml"
STUDY_DIR = DATA_DIR / "studies" / "SOKOLOV2019_HYPOX_MITO"


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def provenance():
    return json.loads((STUDY_DIR / "pride_provenance.json").read_text())


@pytest.fixture(scope="module")
def runs():
    return pd.read_csv(STUDY_DIR / "run_order.tsv", sep="\t")


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_total_proteome_is_three_oysters_per_condition_in_ten_fractions(provenance, runs):
    tp = provenance["total_proteome"]
    assert tp["samples_by_condition"] == {
        "normoxia": ["C2", "C3", "C4"], "hypoxia_24h": ["H1", "H2", "H5"], "reoxygenation_1h": ["R2", "R4", "R5"],
    }
    total = runs[runs.dataset == "total_proteome"]
    assert len(total) == 90 and set(total.groupby("sample_id").size()) == {10}


def test_acquisition_is_blocked_by_condition(provenance):
    order = provenance["total_proteome"]["acquisition_order_of_samples"]
    assert order == ["H1", "C2", "C3", "C4", "R2", "R4", "R5", "H2", "H5"]
    conditions = [s[0] for s in order]
    # Normoxia and reoxygenation are each one consecutive block.
    for c in "CR":
        idx = [i for i, x in enumerate(conditions) if x == c]
        assert idx == list(range(idx[0], idx[0] + 3))


def test_phosphoproteome_has_no_condition_key(provenance, runs):
    pp = provenance["phosphoproteome"]
    assert pp["n_runs"] == 15 and pp["condition_key_deposited"] is False
    assert set(runs.loc[runs.dataset == "phosphoproteome", "condition"]) == {"unknown"}


def test_registered_not_reanalyzed(study, provenance):
    assert study["qc_status"] == "not_started" and study["analysis_status"] == "not_started"
    assert all(c["results_file"] is None for c in study["comparisons"])
    assert study["data_availability"]["processed_data_available"] is False
    assert provenance["pride_files_by_category"] == {"FASTA": 1, "PEAK": 24, "RAW": 105, "RESULT": 24}
    assert not list(STUDY_DIR.glob("*standardized*"))
    for flag in ("low_replication", "batch_confounded", "identifier_mapping_uncertain"):
        assert flag in study["quality_flags"]


def test_comparisons(study):
    comps = {c["comparison_id"]: c for c in study["comparisons"]}
    assert set(comps) == {"hypoxia_24h_vs_normoxia", "reoxygenation_1h_vs_normoxia"}
    assert [cid for cid, c in comps.items() if c.get("meta_analysis_primary")] == ["hypoxia_24h_vs_normoxia"]
    for c in comps.values():
        assert c["stressor_standardized"] == "hypoxia"
        assert (c["sample_size"], c["biological_replicates"]) == (6, 3)


def test_run_order_checksum_matches_provenance(provenance):
    for rel, digest in provenance["outputs_sha256"].items():
        assert hashlib.sha256((REPO_ROOT / rel).read_bytes()).hexdigest() == digest
