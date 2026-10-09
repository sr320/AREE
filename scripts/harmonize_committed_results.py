#!/usr/bin/env python3
"""Harmonize every committed result table of every real study that passed QC.

The Pages workflow and `make real-pool` used to name the real studies to
harmonize, so each newly merged study stayed off the dashboard until someone
edited CI (LUTIER2022_OA_TIPPING was the first to fall through). This walks
the registry instead:

* simulated studies are skipped (the demo studies use the demo crosswalk and
  are harmonized separately);
* studies with `qc_status: failed` are skipped, so their results never reach
  the evidence table (IOCAS2022_OA_ENERGY);
* every comparison with a `results_file` that exists is harmonized, one
  comparison at a time, which works for processed-table and raw_reanalysis
  studies alike. Comparisons with no results_file yet are listed, not errors.

Run with AREE_CROSSWALK pointing at the real crosswalk.

Usage: python scripts/harmonize_committed_results.py [--date YYYY-MM-DD] [--dry-run]
"""
from __future__ import annotations

import argparse
import sys
from datetime import date

from common import REPO_ROOT, STUDIES_DIR, load_yaml


def committed_results() -> tuple[list[tuple[str, str, str]], list[str]]:
    """(study_id, comparison_id, results_file) to harmonize, and notes on what was skipped."""
    jobs, notes = [], []
    for path in sorted(STUDIES_DIR.glob("*.yaml")):
        if path.name.startswith("_"):
            continue
        study = load_yaml(path)
        sid = study["study_id"]
        if study.get("simulated"):
            continue
        if study.get("qc_status") == "failed":
            notes.append(f"{sid}: skipped, qc_status failed")
            continue
        for comp in study.get("comparisons", []):
            results_file = comp.get("results_file")
            if not results_file:
                notes.append(f"{sid}/{comp['comparison_id']}: no results_file yet")
            elif not (REPO_ROOT / results_file).exists():
                notes.append(f"{sid}/{comp['comparison_id']}: results_file missing ({results_file})")
            else:
                jobs.append((sid, comp["comparison_id"], results_file))
    return jobs, notes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--date", default=date.today().isoformat(), help="date_generated (ISO 8601)")
    parser.add_argument("--dry-run", action="store_true", help="list what would be harmonized")
    args = parser.parse_args()

    jobs, notes = committed_results()
    for note in notes:
        print(f"  - {note}")
    if args.dry_run:
        for sid, cid, path in jobs:
            print(f"would harmonize {sid} / {cid} <- {path}")
        return 0

    from harmonize.core import harmonize_processed_table

    for sid, cid, path in jobs:
        df = harmonize_processed_table(sid, REPO_ROOT / path, date_generated=args.date, comparison_id=cid)
        print(f"harmonized {sid} / {cid}: {len(df)} records")
    print(f"{len(jobs)} comparisons harmonized")
    return 0


if __name__ == "__main__":
    sys.exit(main())
