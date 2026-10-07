#!/usr/bin/env python3
"""Check that a raw_reanalysis study's replicates behave like independent animals.

Run after the rnaseq workflow and before harmonizing. Distinct `sample_alias`
values in the deposited metadata do not prove distinct animals: a group whose
"biological replicates" are one RNA sample sequenced several times will pass
every metadata check, then produce DESeq2 output with far too small an lfcSE.
Pooled by inverse variance, such a study would swamp every real one.

For each group, the median per-gene dispersion of size-factor-normalized
counts, (var - mean) / mean^2, is computed over well-expressed genes. Real
biological replicates of M. gigas sit around 0.03-0.05 (CALLA2026_OSHV
Midori control and France groups, 2026-10-01). Values near zero mean
variation is close to Poisson read-sampling noise alone, i.e. technical
replicates. IOCAS2022_OA_ENERGY failed this check at 0.0002-0.0016.

Usage:
  python scripts/check_replicate_dispersion.py \\
      --samplesheet data/studies/STUDY/samplesheet.tsv \\
      --salmon-dir RESULTS/comp_a/rnaseq/salmon [--salmon-dir ...] \\
      [--tx2gene data/reference/GCF_963853765.1/tx2gene.tsv] [--out qc.tsv]

Exits 1 if any group's dispersion is below --min-dispersion.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_TX2GENE = Path("data/reference/GCF_963853765.1/tx2gene.tsv")
# An order of magnitude below the lowest real-replicate value observed so far.
DEFAULT_MIN_DISPERSION = 0.005
MIN_MEAN_COUNT = 500


def gene_counts(salmon_dirs: list[Path], samplesheet: pd.DataFrame, tx2gene: pd.DataFrame) -> pd.DataFrame:
    counts = {}
    for salmon_dir in salmon_dirs:
        for quant in sorted(salmon_dir.glob("*/quant.sf")):
            run = quant.parent.name
            if run not in samplesheet.index:
                continue
            q = pd.read_csv(quant, sep="\t", usecols=["Name", "NumReads"])
            q = q.merge(tx2gene, left_on="Name", right_on="tx")
            counts[samplesheet.loc[run, "sample_id"]] = q.groupby("gene")["NumReads"].sum()
    return pd.DataFrame(counts)


def group_dispersion(counts: pd.DataFrame, groups: pd.Series, min_mean: float = MIN_MEAN_COUNT) -> pd.DataFrame:
    size_factors = counts.sum() / counts.sum().mean()
    norm = counts / size_factors
    rows = []
    for level, members in groups.groupby(groups):
        cols = list(members.index)
        if len(cols) < 2:
            continue
        mean = norm[cols].mean(axis=1)
        var = norm[cols].var(axis=1)
        keep = mean > min_mean
        disp = ((var - mean) / mean**2)[keep]
        rows.append({
            "condition": level,
            "n_samples": len(cols),
            "samples": ",".join(sorted(cols)),
            "n_genes": int(keep.sum()),
            "median_dispersion": round(float(np.median(disp)), 5),
        })
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samplesheet", type=Path, required=True)
    parser.add_argument("--salmon-dir", type=Path, action="append", required=True)
    parser.add_argument("--tx2gene", type=Path, default=DEFAULT_TX2GENE)
    parser.add_argument("--min-dispersion", type=float, default=DEFAULT_MIN_DISPERSION)
    parser.add_argument("--exclude", action="append", default=[], help="sample_id to leave out (repeatable)")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    sheet = pd.read_csv(args.samplesheet, sep="\t", dtype=str).set_index("run_accession")
    tx2gene = pd.read_csv(args.tx2gene, sep="\t", header=None, names=["tx", "gene"], dtype=str)
    counts = gene_counts(args.salmon_dir, sheet, tx2gene)
    counts = counts.drop(columns=[s for s in args.exclude if s in counts.columns])
    if counts.empty:
        print("ERROR: no salmon quant.sf files matched the sample sheet.")
        return 1

    groups = sheet.set_index("sample_id").loc[counts.columns, "condition"]
    table = group_dispersion(counts, groups)
    table["passes"] = table["median_dispersion"] >= args.min_dispersion
    print(table.to_string(index=False))
    if args.out:
        table.to_csv(args.out, sep="\t", index=False)

    if not table["passes"].all():
        print(
            f"\nFAIL: dispersion below {args.min_dispersion} in "
            f"{', '.join(table.loc[~table['passes'], 'condition'])}. "
            "These replicates vary like technical replicates; do not harmonize or pool."
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
