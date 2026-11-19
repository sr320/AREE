"""Compare per-library alt-allele fractions (AF) at shared, well-covered SNPs.

usage: compare_af.py <snps_ad.tsv> <samplesheet.tsv> <out_dir>

Writes pairwise_allele_fractions.tsv and allele_fraction_distribution.tsv.

Expectations:
  distinct individuals   -> genotype classes (hom-ref / het / hom-alt) often
                            disagree; AFs cluster near 0, 0.5 and 1.
  same RNA, re-sequenced -> |dAF| close to read-sampling (binomial) noise,
                            i.e. observed / expected ratio near 1.
  pools of animals       -> many intermediate AFs, not near 0 / 0.5 / 1.

The "expected" column simulates two libraries drawn from one RNA source at
each site's pooled AF and the observed depths (50 replicates, seed 1).
"""
import itertools
import sys

import numpy as np
import pandas as pd

MIN_DP = 30      # depth required in every library
MIN_ALT = 0.15   # alt AF required in at least one library

snps, sheet, out = sys.argv[1], sys.argv[2], sys.argv[3]
raw = pd.read_csv(snps, sep="\t")
raw.columns = [c.split("]")[-1].split(":")[0].strip("#").strip() for c in raw.columns]
samples = [c for c in raw.columns if c.startswith("SRR")]
label = pd.read_csv(sheet, sep="\t").set_index("run_accession")["sample_alias"].to_dict()
samples = sorted(samples, key=lambda s: label[s])

ref = pd.DataFrame({s: raw[s].str.split(",").str[0].replace(".", np.nan).astype(float) for s in samples})
alt = pd.DataFrame({s: raw[s].str.split(",").str[1].replace(".", np.nan).astype(float) for s in samples})
dp = ref + alt
keep = (dp >= MIN_DP).all(axis=1)
af, dp = (alt / dp)[keep], dp[keep]
informative = (af > MIN_ALT).any(axis=1)
af, dp = af[informative], dp[informative]


def gclass(x):
    return np.select([x < 0.1, x > 0.9], [0, 2], 1)


def arm(s):
    return label[s].split(".")[0]


rng = np.random.default_rng(1)
rows = []
for a, b in itertools.combinations(samples, 2):
    pooled = (af[a] * dp[a] + af[b] * dp[b]) / (dp[a] + dp[b])
    na, nb = dp[a].astype(int), dp[b].astype(int)
    expected = np.mean([np.median(np.abs(rng.binomial(na, pooled) / na - rng.binomial(nb, pooled) / nb))
                        for _ in range(50)])
    observed = np.median(np.abs(af[a] - af[b]))
    rows.append({
        "library_a": label[a], "run_a": a, "library_b": label[b], "run_b": b,
        "same_arm": arm(a) == arm(b), "n_sites": len(af),
        "genotype_class_discordance": round((gclass(af[a]) != gclass(af[b])).mean(), 3),
        "median_abs_dAF": round(observed, 4),
        "expected_median_abs_dAF_same_RNA": round(expected, 4),
        "ratio_observed_to_same_RNA": round(observed / expected, 2),
        "AF_correlation": round(np.corrcoef(af[a], af[b])[0, 1], 3),
    })
pd.DataFrame(rows).to_csv(f"{out}/pairwise_allele_fractions.tsv", sep="\t", index=False)

dist = []
for s in samples:
    x = af[s]
    dist.append({
        "library": label[s], "run_accession": s, "median_depth": int(dp[s].median()),
        "af_lt_0_1": round((x < 0.1).mean(), 3),
        "af_0_1_to_0_3": round(((x >= 0.1) & (x < 0.3)).mean(), 3),
        "af_0_3_to_0_7": round(((x >= 0.3) & (x <= 0.7)).mean(), 3),
        "af_0_7_to_0_9": round(((x > 0.7) & (x <= 0.9)).mean(), 3),
        "af_gt_0_9": round((x > 0.9).mean(), 3),
    })
pd.DataFrame(dist).to_csv(f"{out}/allele_fraction_distribution.tsv", sep="\t", index=False)
print(f"{len(af)} sites; wrote pairwise_allele_fractions.tsv and allele_fraction_distribution.tsv to {out}")
