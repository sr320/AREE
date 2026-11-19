"""Print the most expressed transcript of each of the top-N genes by mean TPM."""
import sys
from pathlib import Path

import pandas as pd

sal, tx2gene, n = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3])
tpm = None
for d in sorted(p for p in sal.iterdir() if (p / "quant.sf").exists()):
    q = pd.read_csv(d / "quant.sf", sep="\t", index_col="Name")["TPM"]
    tpm = q if tpm is None else tpm + q
t2g = pd.read_csv(tx2gene, sep="\t", header=None, names=["tx", "gene"]).drop_duplicates("tx").set_index("tx")["gene"]
df = pd.DataFrame({"tpm": tpm})
df["gene"] = t2g.reindex(df.index)
df = df.dropna(subset=["gene"])
best = df.sort_values("tpm", ascending=False).groupby("gene").head(1)
for tx in best.head(n).index:
    print(tx)
