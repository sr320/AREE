"""Sensitivity of the TIMMINS2014_OA_GILL pCO2 comparison (pco2_2800_vs_400).

Two questions, each answered by rerunning workflows/proteomics on the primary
comparison and comparing per-protein log2FC with the committed result:

1. Normalization. The primary run uses NSAF as is (normalization none);
   rerun with --proteomics.normalization median.
2. Run order. Technical replicate 1 of every oyster was acquired in treatment
   blocks (run_order.tsv). Build an NSAF-equivalent matrix from replicates 2
   and 3 only, which were interleaved, and rerun. Per run, spectral counts
   are divided by the run total; protein length cancels in a per-protein
   ratio, so this ranks like NSAF. The authors' per-oyster detection calls are
   kept: a protein absent from an oyster in Table S3 stays absent.

    # 1. write the replicates-2-3 matrix
    python -I sensitivity.py matrix --supplement STAGING/supplementaryFiles.zip --out STAGING/reps23_matrix.tsv
    # 2. run the workflow twice (see the study YAML for the flags), then
    python -I sensitivity.py compare --primary P.tsv --median M.tsv --reps23 R.tsv

`compare` writes sensitivity.tsv and sensitivity_summary.json next to this
script.
"""
from __future__ import annotations

import argparse
import io
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

HERE = Path(__file__).resolve().parent
TABLE_S2 = "12864_2014_7071_MOESM3_ESM.xlsx"


def run_fractions(xlsx: bytes) -> pd.DataFrame:
    """Protein x run spectral-count fractions from Table S2 (one sheet per run)."""
    wb = openpyxl.load_workbook(io.BytesIO(xlsx), read_only=True)
    runs = {}
    for name in wb.sheetnames:
        rows = wb[name].iter_rows(values_only=True)
        header = [str(h) for h in next(rows)]
        i_prot, i_spc = header.index("protein"), header.index("tot indep spectra")
        counts = {}
        for r in rows:
            if r[i_prot] is None:
                continue
            for prot in str(r[i_prot]).split(","):
                counts[prot] = max(counts.get(prot, 0), float(r[i_spc] or 0))
        s = pd.Series(counts)
        runs[name] = s / s.sum()
    return pd.DataFrame(runs).fillna(0)


def build_matrix(supplement: Path, out: Path) -> None:
    with zipfile.ZipFile(supplement) as zf:
        fractions = run_fractions(zf.read(TABLE_S2))
    nsaf = pd.read_csv(HERE / "nsaf_matrix.tsv", sep="\t").set_index("protein_accession")
    reps23 = pd.DataFrame({
        col: fractions[[f"{int(col[2:])}_02", f"{int(col[2:])}_03"]].mean(axis=1).reindex(nsaf.index).fillna(0)
        for col in nsaf.columns
    })
    reps23.where(nsaf > 0, 0).to_csv(out, sep="\t")
    print(f"wrote {out}")


def compare(primary: Path, median: Path, reps23: Path) -> None:
    def load(p):
        return pd.read_csv(p, sep="\t").set_index("protein_accession")[["log2FC", "lfcSE", "padj"]]

    base = load(primary)
    table = base.rename(columns=lambda c: f"{c}_primary")
    summary = {}
    for label, path in (("median_normalization", median), ("replicates_2_3", reps23)):
        other = load(path)
        table = table.join(other.rename(columns=lambda c, lab=label: f"{c}_{lab}"), how="outer")
        both = base.join(other, rsuffix="_o", how="inner")
        summary[label] = {
            "n_proteins_tested_in_both": int(len(both)),
            "log2FC_pearson_r": round(float(np.corrcoef(both["log2FC"], both["log2FC_o"])[0, 1]), 3),
            "share_same_direction": round(float((np.sign(both["log2FC"]) == np.sign(both["log2FC_o"])).mean()), 3),
            "mean_log2FC_shift": round(float((both["log2FC_o"] - both["log2FC"]).mean()), 3),
            "n_padj_below_0.05": int((other["padj"] < 0.05).sum()),
        }
    table.sort_index().to_csv(HERE / "sensitivity.tsv", sep="\t", float_format="%.6g")
    (HERE / "sensitivity_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("matrix")
    m.add_argument("--supplement", type=Path, required=True)
    m.add_argument("--out", type=Path, required=True)
    c = sub.add_parser("compare")
    c.add_argument("--primary", type=Path, required=True)
    c.add_argument("--median", type=Path, required=True)
    c.add_argument("--reps23", type=Path, required=True)
    args = ap.parse_args()
    if args.cmd == "matrix":
        build_matrix(args.supplement, args.out)
    else:
        compare(args.primary, args.median, args.reps23)
    return 0


if __name__ == "__main__":
    sys.exit(main())
