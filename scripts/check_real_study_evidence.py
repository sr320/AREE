#!/usr/bin/env python3
"""Assert that harmonizing the real study against the real crosswalk still works.

Run after `aree harmonize --study HESSER2024_VCOR`. Guards the properties that
distinguish the real-data path from the demo path, and that would otherwise
regress silently:

* the expected number of evidence records is produced;
* nothing from a real study is ever flagged `simulated`;
* identifier resolution against the real NCBI/UniProt crosswalk has not decayed;
* the study still reports no standard error, so nothing downstream has begun
  imputing one (see docs/first_real_study.md).

Usage: python scripts/check_real_study_evidence.py [--study STUDY_ID]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

EVIDENCE_TABLE = Path("reports/evidence/evidence_table.tsv")

EXPECTED_RECORDS = 351
MIN_RESOLVED_FRACTION = 0.80


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", default="HESSER2024_VCOR")
    args = parser.parse_args()

    if not EVIDENCE_TABLE.exists():
        print(f"ERROR: no evidence table at {EVIDENCE_TABLE}; run `aree harmonize` first.")
        return 1

    df = pd.read_csv(EVIDENCE_TABLE, sep="\t")
    real = df[df["study_id"] == args.study]

    failures = []
    if len(real) != EXPECTED_RECORDS:
        failures.append(f"expected {EXPECTED_RECORDS} evidence records, got {len(real)}")

    if real["simulated"].astype(str).str.lower().eq("true").any():
        failures.append("a real study is flagged `simulated`")

    resolved = (real["mapping_confidence"] != "unresolved").mean()
    if resolved < MIN_RESOLVED_FRACTION:
        failures.append(
            f"identifier resolution fell to {resolved:.1%} (floor {MIN_RESOLVED_FRACTION:.0%})"
        )

    # The published source reports neither a standard error nor an unadjusted
    # p-value. If either becomes non-empty, something is imputing statistics the
    # source never reported — which AREE must never do.
    for column in ("standard_error", "p_value"):
        if real[column].notna().any():
            failures.append(
                f"`{column}` is populated for {args.study}, but the source reports none — "
                "a statistic is being imputed"
            )

    if failures:
        print(f"FAIL: {args.study}")
        for f in failures:
            print(f"  - {f}")
        return 1

    print(
        f"OK: {args.study} — {len(real)} evidence records, "
        f"{resolved:.1%} identifiers resolved, no imputed statistics"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
