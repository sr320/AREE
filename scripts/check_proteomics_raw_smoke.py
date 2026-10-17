"""Check the proteomics raw-mode CI run on tests/fixtures/proteomics_raw.

The fixture plants two effects (PX_up +2.0, PX_down -1.5), one protein
detected only in treatment and one never detected, in a linear NSAF-like
matrix with two extra columns outside the sample sheet. The run must recover
the effects in direction and rough size, give every tested row a standard
error and a spectral-count analysis_method, keep untestable proteins out of the
results, compare 4 vs 4 samples, and record its options in the manifest.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

PLANTED = {"PX_up": 2.0, "PX_down": -1.5}
NOT_TESTED = {"PX_treatonly": "detected_too_rarely_in_control", "PX_absent": "not_detected"}


def main(outdir: str) -> int:
    root = Path(outdir) / "proteomics"
    stem = "CI_PROT_treatment_vs_control"
    res = pd.read_csv(root / "standardize_output" / f"{stem}_protein_abundance_standardized.tsv", sep="\t")
    status = pd.read_csv(root / "differential_abundance" / f"{stem}_detection_status.tsv", sep="\t")
    manifest = json.loads((root / "manifest" / f"{stem}_manifest.json").read_text())
    errors = []

    by_id = res.set_index("protein_accession")
    for pid, lfc in PLANTED.items():
        row = by_id.loc[pid]
        if not (row.padj < 0.05 and abs(row.log2FC - lfc) < 0.5):
            errors.append(f"{pid}: log2FC {row.log2FC:.2f}, padj {row.padj:.2g}; planted {lfc:+.1f}")
    if res["lfcSE"].isna().any() or (res["lfcSE"] <= 0).any():
        errors.append("every tested row needs a positive lfcSE")
    if set(res["analysis_method"]) != {"spectral_count_limma_moderated_t"}:
        errors.append(f"analysis_method: {sorted(set(res['analysis_method']))}")
    leaked = set(NOT_TESTED) & set(res["protein_accession"])
    if leaked:
        errors.append(f"untestable proteins in the results: {sorted(leaked)}")

    st = status.set_index("protein_accession")
    for pid, expected in NOT_TESTED.items():
        if st["status"].get(pid) != expected:
            errors.append(f"detection status {pid}: {st['status'].get(pid)!r}, expected {expected!r}")
    if set(st["n_replicates_control"]) != {4} or set(st["n_replicates_treatment"]) != {4}:
        errors.append("the sample sheet names 4 vs 4 samples; other_1/other_2 must be left out")

    params = manifest["parameters"]["proteomics"]
    for key, value in {"input_scale": "linear", "normalization": "none", "quantification": "spectral_count"}.items():
        if params.get(key) != value:
            errors.append(f"manifest parameters.proteomics.{key} = {params.get(key)!r}")
    da = manifest["qc_metrics"]["differential_abundance_qc"]
    if not da or da["n_proteins_tested"] != len(res):
        errors.append("manifest qc_metrics.differential_abundance_qc does not match the results")
    if len(manifest["inputs"]["input_files_sha256"]) != 2:
        errors.append("manifest should checksum the abundance matrix and the sample sheet")

    for e in errors:
        print(f"FAIL: {e}")
    print(f"{len(res)} tested, {int((st['status'] != 'tested').sum())} not tested, {len(errors)} problem(s)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
