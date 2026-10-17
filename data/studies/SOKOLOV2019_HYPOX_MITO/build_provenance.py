"""Record what PXD011365 holds, and why it is registered but not reanalyzed.

PRIDE PXD011365 (CC0) has two datasets from gill mitochondria of oysters under
normoxia (C), 24 h hypoxia (H) and 1 h reoxygenation after hypoxia (R):

- total mitochondrial proteome: 9 oysters (C2-C4, H1/H2/H5, R2/R4/R5), each run
  as 10 GeLC-MS fractions on an LTQ Orbitrap Classic, with one Scaffold mzid
  export per oyster ("Mudpit_..._1.mzid.gz");
- enriched phosphoproteome: 15 runs (PP_01-PP_15) on an LTQ Orbitrap Velos,
  5 per condition according to the PRIDE protocol, with no key saying which
  run is which condition.

    python -I data/studies/SOKOLOV2019_HYPOX_MITO/build_provenance.py --download STAGING_DIR

fetches the PRIDE project and file records and the first KB of every .raw
file, reads each run's acquisition start (Thermo header, FILETIME at byte 40)
and writes

  data/studies/SOKOLOV2019_HYPOX_MITO/run_order.tsv
  data/studies/SOKOLOV2019_HYPOX_MITO/pride_provenance.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import struct
import sys
import urllib.request
from collections import Counter
from datetime import date, datetime, timedelta
from pathlib import Path

import pandas as pd

ACCESSION = "PXD011365"
STUDY_ID = "SOKOLOV2019_HYPOX_MITO"
PRIDE_API = f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{ACCESSION}"
PRIDE_FTP = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/11/PXD011365"

REPO_ROOT = Path(__file__).resolve().parents[3]
OUT = REPO_ROOT / "data" / "studies" / STUDY_ID
CONDITION = {"C": "normoxia", "H": "hypoxia_24h", "R": "reoxygenation_1h"}


def fetch(url: str, byte_range: str | None = None) -> bytes:
    headers = {"User-Agent": "AREE-curation"}
    if byte_range:
        headers["Range"] = f"bytes={byte_range}"
    with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=300) as resp:
        return resp.read()


def raw_start_time(header: bytes) -> str:
    ticks = struct.unpack_from("<Q", header, 40)[0]
    return (datetime(1601, 1, 1) + timedelta(microseconds=ticks / 10)).isoformat(timespec="seconds")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--download", type=Path, required=True, help="staging directory for fetched files")
    args = ap.parse_args()
    args.download.mkdir(parents=True, exist_ok=True)

    project_raw, files_raw = fetch(PRIDE_API), fetch(f"{PRIDE_API}/files/all")
    (args.download / "pride_project.json").write_bytes(project_raw)
    (args.download / "pride_files.json").write_bytes(files_raw)
    project, files = json.loads(project_raw), json.loads(files_raw)

    runs = []
    for f in sorted(x["fileName"] for x in files if x["fileCategory"]["value"] == "RAW"):
        stem = f.removesuffix(".raw")
        m = re.search(r"_StM_([CHR])(\d)_(\d+)$", stem)
        if m:
            dataset, sample, condition, fraction = ("total_proteome", f"{m.group(1)}{m.group(2)}",
                                                    CONDITION[m.group(1)], int(m.group(3)))
        else:
            m = re.search(r"_PP_(\d+)$", stem)
            if not m:
                sys.exit(f"unrecognised run name {f}")
            dataset, sample, condition, fraction = "phosphoproteome", f"PP_{m.group(1)}", "unknown", 1
        runs.append({"raw_file": f, "dataset": dataset, "sample_id": sample, "condition": condition,
                     "fraction": fraction, "acquired": raw_start_time(fetch(f"{PRIDE_FTP}/{f}", "0-1023"))})
    runs.sort(key=lambda r: r["acquired"])
    for i, r in enumerate(runs, start=1):
        r["run_order"] = i
    table = pd.DataFrame(runs)
    table.to_csv(OUT / "run_order.tsv", sep="\t", index=False)

    total = table[table.dataset == "total_proteome"]
    sample_order = list(dict.fromkeys(total["sample_id"]))
    if total.groupby("sample_id")["fraction"].nunique().ne(10).any():
        sys.exit("expected 10 fractions per total-proteome sample")
    # Each sample's fractions must be contiguous in the run order (one sample at a time).
    for s, g in total.groupby("sample_id"):
        if g["run_order"].max() - g["run_order"].min() != 9:
            sys.exit(f"sample {s}: fractions are not acquired consecutively")
    phospho_order = list(table.loc[table.dataset == "phosphoproteome", "sample_id"])

    provenance = {
        "study_id": STUDY_ID,
        "accession": ACCESSION,
        "repository": "PRIDE",
        "license": project["license"],
        "submission_type": project["submissionType"],
        "date_retrieved": date.today().isoformat(),
        "generated_by": f"data/studies/{STUDY_ID}/build_provenance.py",
        "source_urls": [PRIDE_API, f"{PRIDE_API}/files/all", f"{PRIDE_FTP}/<run>.raw (bytes 0-1023)"],
        "pride_reference": project["references"][0]["referenceLine"],
        "pride_instrument": [i["name"] for i in project["instruments"]],
        "pride_files_by_category": dict(sorted(Counter(f["fileCategory"]["value"] for f in files).items())),
        "total_proteome": {
            "n_samples": len(sample_order),
            "samples_by_condition": {c: sorted(total.loc[total.condition == c, "sample_id"].unique())
                                     for c in CONDITION.values()},
            "fractions_per_sample": 10,
            "acquisition_order_of_samples": sample_order,
            "acquisition_window": [total["acquired"].min(), total["acquired"].max()],
            "quantitative_files": "one Scaffold 4.8.4 mzIdentML export per sample (Mudpit_*.mzid.gz); "
                                  "no per-sample abundance table deposited",
        },
        "phosphoproteome": {
            "n_runs": len(phospho_order),
            "acquisition_order": phospho_order,
            "acquisition_window": [table.loc[table.dataset == "phosphoproteome", "acquired"].min(),
                                   table.loc[table.dataset == "phosphoproteome", "acquired"].max()],
            "condition_key_deposited": False,
        },
        "outputs_sha256": {f"data/studies/{STUDY_ID}/run_order.tsv":
                           hashlib.sha256((OUT / "run_order.tsv").read_bytes()).hexdigest()},
    }
    (OUT / "pride_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"wrote run_order.tsv ({len(runs)} runs) and pride_provenance.json; sample order {sample_order}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
