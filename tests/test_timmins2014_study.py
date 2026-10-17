"""Curation tests for TIMMINS2014_OA_GILL.

Real proteomics study (PRIDE PXD000835; Timmins-Schiffman et al. 2014, BMC
Genomics), reanalyzed from the paper's per-oyster NSAF table with
workflows/proteomics. These tests pin the design (4 individual oysters per
pCO2 x mechanical-stress group, 3 technical runs each), the data-checked
column -> oyster mapping, the run order, the reanalysis settings, the
low-power null result, and the sensitivity checks the QC decision rests on.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter

import pandas as pd
import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_PATH = STUDIES_DIR / "TIMMINS2014_OA_GILL.yaml"
STUDY_DIR = DATA_DIR / "studies" / "TIMMINS2014_OA_GILL"
COMPARISONS = ("pco2_2800_vs_400", "mech_stress_vs_none_at_400", "mech_stress_vs_none_at_2800")
GROUPS = {("2800", "no"): {2, 5, 8, 11}, ("2800", "yes"): {26, 29, 32, 35},
          ("400", "no"): {221, 224, 227, 230}, ("400", "yes"): {242, 245, 248, 251}}


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def provenance():
    return json.loads((STUDY_DIR / "inputs_provenance.json").read_text())


@pytest.fixture(scope="module")
def samples():
    return pd.read_csv(STUDY_DIR / "samplesheet.tsv", sep="\t", dtype=str)


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_sixteen_individual_oysters_in_four_groups(samples, provenance):
    assert len(samples) == provenance["n_oysters"] == 16
    got = {(p, m): set(g["oyster_number"].astype(int)) for (p, m), g in samples.groupby(["pco2_uatm", "mechanical_stress"])}
    assert got == GROUPS


def test_nsaf_columns_are_mapped_to_oysters_by_data(samples, provenance):
    # Each Table S3 column matches its own oyster's spectral counts far better
    # than any other oyster's; a column-order mistake would break this.
    r = samples["pearson_check_r"].astype(float)
    nxt = samples["next_best_r"].astype(float)
    assert (r > 0.5).all() and (nxt < 0.2).all()
    assert provenance["min_mapping_margin_r"] > 0.4
    nsaf = pd.read_csv(STUDY_DIR / "nsaf_matrix.tsv", sep="\t")
    assert list(nsaf.columns[1:]) == list(samples["sample_id"])
    assert len(nsaf) == provenance["n_proteins"] == 1616
    assert nsaf["protein_accession"].str.fullmatch(r"CGI_\d+").all()


def test_comparison_sheets_are_four_vs_four(samples):
    group_of = {(row.pco2_uatm, row.mechanical_stress): set() for row in samples.itertuples()}
    for row in samples.itertuples():
        group_of[(row.pco2_uatm, row.mechanical_stress)].add(row.sample_id)
    expected = {
        "pco2_2800_vs_400": (("400", "no"), ("2800", "no")),
        "mech_stress_vs_none_at_400": (("400", "no"), ("400", "yes")),
        "mech_stress_vs_none_at_2800": (("2800", "no"), ("2800", "yes")),
    }
    for cid, (ctl, trt) in expected.items():
        sheet = pd.read_csv(STUDY_DIR / f"samplesheet_{cid}.tsv", sep="\t")
        assert set(sheet.loc[sheet.group == "control", "sample_id"]) == group_of[ctl]
        assert set(sheet.loc[sheet.group == "treatment", "sample_id"]) == group_of[trt]


def test_run_order_blocks_only_technical_replicate_1(provenance):
    runs = pd.read_csv(STUDY_DIR / "run_order.tsv", sep="\t", dtype=str)
    assert len(runs) == 48 and runs["run_order"].astype(int).tolist() == list(range(1, 49))
    assert Counter(runs["sample_id"]).most_common()[-1][1] == 3
    assert provenance["technical_replicate_1_order"] == ["400"] * 8 + ["2800"] * 8
    # "101B"/"103B" in the file names is the pCO2 split, not an MS batch.
    assert set(zip(runs["sample_label"], runs["pco2_uatm"])) == {("101B", "2800"), ("103B", "400")}


def test_reanalysis_settings_in_every_manifest():
    for cid in COMPARISONS:
        m = json.loads((STUDY_DIR / f"TIMMINS2014_OA_GILL_{cid}_manifest.json").read_text())
        p = m["parameters"]["proteomics"]
        assert m["mode"] == "raw_reanalysis" and p["comparison_id"] == cid
        assert (p["input_scale"], p["normalization"], p["quantification"], int(p["min_detected_per_arm"])) == \
            ("linear", "none", "spectral_count", 3)
        da = m["qc_metrics"]["differential_abundance_qc"]
        assert (da["n_replicates_control"], da["n_replicates_treatment"]) == (4, 4)
        assert da["n_proteins_in"] == 1616


def test_results_are_a_low_power_null_with_standard_errors(study):
    for comp in study["comparisons"]:
        res = pd.read_csv(REPO_ROOT / comp["results_file"], sep="\t")
        assert 600 <= len(res) <= 650
        assert res["lfcSE"].notna().all() and (res["lfcSE"] > 0).all()
        assert set(res["analysis_method"]) == {"spectral_count_limma_moderated_t"}
        assert (res["padj"] < 0.05).sum() == 0
        # Close to the null: a few percent of p-values below 0.05.
        assert (res["pvalue"] < 0.05).mean() < 0.1


def test_every_protein_has_a_detection_status():
    for cid in COMPARISONS:
        st = pd.read_csv(STUDY_DIR / f"TIMMINS2014_OA_GILL_{cid}_detection_status.tsv", sep="\t")
        res = pd.read_csv(STUDY_DIR / f"TIMMINS2014_OA_GILL_{cid}_protein_abundance_standardized.tsv", sep="\t")
        assert len(st) == 1616
        assert set(st.loc[st.status == "tested", "protein_accession"]) == set(res["protein_accession"])


def test_sensitivity_supports_the_qc_decision(study):
    s = json.loads((STUDY_DIR / "sensitivity_summary.json").read_text())
    assert s["median_normalization"]["log2FC_pearson_r"] > 0.95
    assert s["replicates_2_3"]["log2FC_pearson_r"] > 0.8
    assert abs(s["replicates_2_3"]["mean_log2FC_shift"]) < 0.1
    assert study["qc_status"] == "passed_with_warnings"


def test_primary_comparison_and_classifications(study):
    comps = {c["comparison_id"]: c for c in study["comparisons"]}
    assert [cid for cid, c in comps.items() if c.get("meta_analysis_primary")] == ["pco2_2800_vs_400"]
    assert comps["pco2_2800_vs_400"]["stressor_standardized"] == "ocean_acidification"
    assert comps["pco2_2800_vs_400"]["resilience_classification"] == "resilience"
    for cid in ("mech_stress_vs_none_at_400", "mech_stress_vs_none_at_2800"):
        assert comps[cid]["resilience_classification"] == "exposure_only"
        assert comps[cid]["stressor_standardized"] == "handling_transport_stress"
    for c in comps.values():
        assert (c["sample_size"], c["biological_replicates"]) == (8, 4)


def test_outputs_match_provenance(provenance):
    for rel, digest in provenance["outputs_sha256"].items():
        assert hashlib.sha256((REPO_ROOT / rel).read_bytes()).hexdigest() == digest, rel
