"""Curation tests for ZHONG2025_PONV1_MORTALITY.

Adult C. gigas sampled in Baynes Sound, BC, during a 2020 mass-mortality event
(Zhong et al. 2025): 96 individually sequenced mantle and hemolymph libraries
from diseased farmed and apparently healthy wild oysters. These tests pin the
decisions that are easy to erode later: farm vs wild is registered only as a
confounded exposure contrast, never as pooling-grade evidence; the within-farm
PONV1-positive vs PONV1-negative comparison stays unregistered, because the
rule fixed before any viral counts gave fewer than 5 positive animals; the
PONV1 calls rest on both genome segments; and the sub-threshold PONV1 in wild
E4 and D5 is recorded as a conflict with the paper rather than dropped.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter

import pytest
import yaml

from common import DATA_DIR, REPO_ROOT, STUDIES_DIR

STUDY_ID = "ZHONG2025_PONV1_MORTALITY"
STUDY_PATH = STUDIES_DIR / f"{STUDY_ID}.yaml"
STUDY_DIR = DATA_DIR / "studies" / STUDY_ID
SEGMENT_DIR = STUDY_DIR / "ponv1_per_segment"
PONV1_DIR = STUDY_DIR / "pathogen_reference"
PONV1_FASTA = "PONV1_PQ030832.1_PV579160.1.fasta"
FARM = {"StnI", "StnL", "StnM"}
COMPARISONS = {
    "farm_mortality_vs_wild_mantle": ("Mantle", 27, 26),
    "farm_mortality_vs_wild_hemolymph": ("Hemolymph", 19, 24),
}
# rna.fna and genomic.gtf of the provenance-checked GCF_963853765.1 bundle
# (as recorded by DUPERRET2025_HEAT_RNASEQ's manifests).
BUNDLE_SHA256 = {
    "rna.fna": "297d49baac2fd3187db957c3dcd150e58ed2a61c3867b61dbdd83ba9b683c386",
    "genomic.gtf": "ab5443f889932ef3e860d47e23801dd94a2b6a4a508f807e0349dfaccea79600",
}


@pytest.fixture(scope="module")
def study():
    return yaml.safe_load(STUDY_PATH.read_text())


@pytest.fixture(scope="module")
def comparisons(study):
    return {c["comparison_id"]: c for c in study["comparisons"]}


@pytest.fixture(scope="module")
def samplesheet():
    return _read_tsv(STUDY_DIR / "samplesheet.tsv")


def _read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t"))


def _manifest(cid):
    return json.loads((STUDY_DIR / f"{STUDY_ID}_{cid}_workflow_manifest.json").read_text())


def _segments(tissue):
    return _read_tsv(SEGMENT_DIR / f"ponv1_per_segment_{tissue}.tsv")


def test_study_validates():
    from intake.schema_validate import validate_study_file

    result = validate_study_file(STUDY_PATH)
    assert result.valid, result.errors


def test_one_biosample_per_animal_per_tissue(samplesheet):
    assert len(samplesheet) == 96
    assert len({r["sample_accession"] for r in samplesheet}) == 96
    arms = Counter(("farm" if r["station"] in FARM else "wild", r["tissue_label"]) for r in samplesheet)
    assert arms == {("farm", "Mantle"): 27, ("wild", "Mantle"): 26,
                    ("farm", "Hemolymph"): 19, ("wild", "Hemolymph"): 24}


def test_farm_vs_wild_is_a_confounded_exposure_contrast(comparisons):
    """Farm vs wild differs in genetics, husbandry and site as well as disease:
    recorded, never flagged for pooling."""
    assert set(comparisons) == set(COMPARISONS)
    for c in comparisons.values():
        assert not c.get("meta_analysis_primary")
        assert c["contrast_type"] == "exposure_response"
        assert c["stressor_standardized"] == "multi_stressor"
        assert c["phenotype"] == "disease_susceptibility"
        assert c["resilience_classification"] == "disease"
        assert c["life_stage"] == "adult"
        assert "Etiology not established" in c["stressor_original"]


def test_no_within_farm_ponv1_comparison_is_registered(comparisons):
    """The rule fixed on 2026-10-14 needs at least 5 PONV1-positive farm mantle
    animals; there are 4."""
    assert not any("ponv1" in cid.lower() for cid in comparisons)
    calls = Counter((r["arm"], r["prespecified_call"]) for r in _segments("mantle"))
    assert calls["farm", "positive"] == 4 < 5


@pytest.mark.parametrize("tissue,farm_positive,farm_n,wild_n", [
    ("mantle", {"L10", "I8", "I16", "I6"}, 27, 26),
    ("hemolymph", {"L10"}, 19, 24),
])
def test_prespecified_ponv1_calls(tissue, farm_positive, farm_n, wild_n):
    rows = _segments(tissue)
    assert Counter(r["arm"] for r in rows) == {"farm": farm_n, "wild": wild_n}
    positive = {r["animal"] for r in rows if r["prespecified_call"] == "positive"}
    assert positive == farm_positive
    for r in rows:
        meets = float(r["ponv1_per_million"]) >= 100 and r["both_segments_detected"] == "True"
        assert (r["prespecified_call"] == "positive") == meets
        both = int(r["seg1_PQ030832_pairs"]) > 0 and int(r["seg2_PV579160_pairs"]) > 0
        assert (r["both_segments_detected"] == "True") == both


def test_wild_ponv1_signal_is_recorded_not_explained_by_index_hopping(study):
    """Wild E4 carries PONV1 in both tissues, sequenced on different lanes;
    lane-3 wild hemolymph beside 16,640 farm PONV1 pairs shows at most 3."""
    mantle = {r["animal"]: r for r in _segments("mantle") if r["arm"] == "wild"}
    hemo = {r["animal"]: r for r in _segments("hemolymph") if r["arm"] == "wild"}
    assert int(mantle["E4"]["ponv1_pairs"]) == 2575
    assert int(mantle["D5"]["ponv1_pairs"]) == 740
    assert int(hemo["E4"]["ponv1_pairs"]) == 188
    assert mantle["E4"]["flowcell_lane"] != hemo["E4"]["flowcell_lane"]
    lane3 = [r for r in _segments("hemolymph") if r["flowcell_lane"] == "HG5LKDSX2:3"]
    farm_pairs = sum(int(r["ponv1_pairs"]) for r in lane3 if r["arm"] == "farm")
    other_wild = [int(r["ponv1_pairs"]) for r in lane3 if r["arm"] == "wild" and r["animal"] != "E4"]
    assert farm_pairs == 16_640
    assert len(other_wild) == 23 and max(other_wild) <= 3
    notes = study["provenance"]["curation_notes"]
    assert "Index" in notes and "hopping therefore does not explain the wild signal" in notes
    assert "E4" in study["limitations"]


def test_ponv1_reference_matches_its_provenance():
    prov = json.loads((PONV1_DIR / "PONV1_PQ030832.1_PV579160.1.provenance.json").read_text())
    data = (PONV1_DIR / PONV1_FASTA).read_bytes()
    assert prov["artifact"] == PONV1_FASTA
    assert hashlib.sha256(data).hexdigest() == prov["artifact_sha256"]
    assert prov["accession"] == ["PQ030832.1", "PV579160.1"]
    records = [r.split("\n", 1) for r in data.decode().split(">")[1:]]
    assert [header.split()[0] for header, _ in records] == ["PQ030832.1", "PV579160.1"]
    assert [len(seq.replace("\n", "")) for _, seq in records] == [41_709, 22_622]


@pytest.mark.parametrize("cid", COMPARISONS)
def test_manifest_records_the_run(cid):
    tissue, n_farm, n_wild = COMPARISONS[cid]
    manifest = _manifest(cid)
    params = manifest["parameters"]
    assert params["viral_reference"].endswith(PONV1_FASTA)
    assert params["exclude_genes"].endswith("GCF_963853765.1-RS_2024_06_rRNA.txt")
    assert params["exclude_samples"] == []
    assert all(level.endswith(f"_{tissue}") for level in params["control_levels"] + params["treatment_levels"])
    assert {lvl.split("_")[0] for lvl in params["treatment_levels"]} <= FARM
    assert not {lvl.split("_")[0] for lvl in params["control_levels"]} & FARM
    design = manifest["deseq2_design"]
    assert design["n_libraries_or_units"] == n_farm + n_wild
    assert design["residual_df"] == n_farm + n_wild - 2
    assert design["exclude_genes"]["n_listed"] == 68
    files = {f["name"]: f["sha256"] for f in manifest["reference"]["files"]}
    for name, sha in BUNDLE_SHA256.items():
        assert files[name] == sha
    assert manifest["software_versions"]["SALMON_QUANT"]["salmon"] == "2.6.0"
    assert manifest["software_versions"]["TRIM_FASTP"]["fastp"] == "1.3.6"


@pytest.mark.parametrize("cid,expected", [
    ("farm_mortality_vs_wild_mantle", (31_813, 17_647, 13_323)),
    ("farm_mortality_vs_wild_hemolymph", (31_807, 18_768, 17_150)),
])
def test_results_carry_standard_errors(comparisons, cid, expected):
    rows = _read_tsv(REPO_ROOT / comparisons[cid]["results_file"])
    n, n_sig, n_up = expected
    assert len(rows) == n
    assert all(r["lfcSE"] not in ("", "NA") for r in rows)
    sig = [r for r in rows if r["padj"] not in ("", "NA") and float(r["padj"]) < 0.05]
    assert len(sig) == n_sig
    assert sum(float(r["log2FoldChange"]) > 0 for r in sig) == n_up


def test_caveats_stay_on_record(study):
    assert {"unbalanced_design", "single_timepoint", "raw_data_available"} <= set(study["quality_flags"])
    limitations = study["limitations"]
    for phrase in ("confounded with disease state", "Etiology is not established",
                   "16-36%", "partly", "without rRNA depletion"):
        assert phrase in limitations, phrase
    assert study["qc_status"] == "passed_with_warnings"
    assert study["analysis_status"] == "harmonized"
    assert "CONFIRMED" in STUDY_PATH.read_text()
