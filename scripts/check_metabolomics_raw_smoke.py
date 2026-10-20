"""Check the metabolomics raw-mode CI run on tests/fixtures/metabolomics_raw.

The fixture plants three effects (LC_up +2.0, LC_down -1.5, GC_up +1.0), one
metabolite detected only in treatment and one never detected. The run must
recover the effects in direction and rough size, give every tested row a
standard error, keep internal standards and undetected metabolites out of the
results, and list every non-standard feature in the missingness report.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

PLANTED = {"LC_up": 2.0, "LC_down": -1.5, "GC_up": 1.0}
NOT_TESTED = {"LC_treatonly": "detected_too_rarely_in_control", "LC_lod": "not_detected"}
INTERNAL_STANDARDS = {"LC_IS", "GC_IS"}


def main(outdir: str) -> int:
    root = Path(outdir) / "metabolomics"
    stem = "CI_METAB_treatment_vs_control"
    res = pd.read_csv(root / "standardized" / f"{stem}_metabolite_features_standardized.tsv", sep="\t")
    miss = pd.read_csv(root / "differential_abundance" / f"{stem}_missingness.tsv", sep="\t")
    manifest = json.loads((root / "manifest" / f"{stem}_manifest.json").read_text())
    errors = []

    by_id = res.set_index("feature_id")
    for fid, lfc in PLANTED.items():
        row = by_id.loc[fid]
        if not (row.padj < 0.05 and abs(row.log2FC - lfc) < 0.5):
            errors.append(f"{fid}: log2FC {row.log2FC:.2f}, padj {row.padj:.2g}; planted {lfc:+.1f}")
    if res["lfcSE"].isna().any() or (res["lfcSE"] <= 0).any():
        errors.append("every tested row needs a positive lfcSE")
    if set(res["analysis_method"]) != {"targeted_MS_limma_moderated_t"}:
        errors.append(f"analysis_method: {sorted(set(res['analysis_method']))}")
    leaked = (set(NOT_TESTED) | INTERNAL_STANDARDS) & set(res["feature_id"])
    if leaked:
        errors.append(f"untestable features in the results: {sorted(leaked)}")

    status = miss.set_index("feature_id")["status"]
    for fid, expected in NOT_TESTED.items():
        if status.get(fid) != expected:
            errors.append(f"missingness {fid}: {status.get(fid)!r}, expected {expected!r}")
    if INTERNAL_STANDARDS & set(miss["feature_id"]):
        errors.append("internal standards belong in normalization QC, not the missingness report")
    if set(miss["n_replicates_control"]) != {3}:
        errors.append("replicate_unit tank should give 3 control units")

    params = manifest["parameters"]["metabolomics"]
    for key, value in {"normalization": "internal_standard", "replicate_unit": "tank"}.items():
        if params.get(key) != value:
            errors.append(f"manifest parameters.metabolomics.{key} = {params.get(key)!r}")

    for e in errors:
        print(f"FAIL: {e}")
    print(f"{len(res)} tested, {int((status != 'tested').sum())} not tested, {len(errors)} problem(s)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
