"""Build the sample sheet and provenance for CRANDALL2022_HEAT_SEED (PXD015434).

PRIDE PXD015434 (CC0) deposits 16 DIA runs (.raw) and their EncyclopeDIA
.elib files, with no sample metadata. The authors' Zenodo release
(10.5281/zenodo.5708415, grace-ac/paper-pacific.oyster-larvae v2.1.0) holds
what ties runs to samples:

  data/20190403-MSstats-report.csv.zip   Skyline export: Condition and
                                         BioReplicate for every run file
  data/20190403-2015Cgseed-protcomp.csv  the MSstats groupComparison result
  data/2015X15P8LargescaleSeed29CtoAmbient.xlsx
                                         the hatchery notebook, including the
                                         proteomic sample key (Sheet1)

The Skyline replicate names (Oysterseed_01, _02, _13, _14; Zenodo
images/skyline-images/40-condition-and-newbioreps.PNG) are hatchery sample
numbers 1, 2, 13 and 14 in that key: the 23 C control silo (Silo 9) and the
29 C control silo (Silo 8) on 11 and 14 Sep 2015, exposure days 5 and 8. That
link comes from a screenshot, so it is written down here as SKYLINE_REPLICATES
and checked against the key and the MSstats report rather than parsed.

    python -I data/studies/CRANDALL2022_HEAT_SEED/build_samplesheet.py --download STAGING_DIR

fetches the PRIDE records and the Zenodo zip (about 68 MB, md5-checked against
the Zenodo API) into STAGING_DIR and writes

  data/studies/CRANDALL2022_HEAT_SEED/samplesheet.tsv      one row per DIA run
  data/studies/CRANDALL2022_HEAT_SEED/pride_provenance.json

Nothing from the Zenodo release is committed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sys
import urllib.request
import zipfile
from collections import Counter
from datetime import date, datetime
from pathlib import Path

import openpyxl

ACCESSION = "PXD015434"
STUDY_ID = "CRANDALL2022_HEAT_SEED"
PRIDE_API = f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{ACCESSION}"
ZENODO_RECORD = 5708415
ZENODO_API = f"https://zenodo.org/api/records/{ZENODO_RECORD}"
ZIP_ROOT = "grace-ac-paper-pacific.oyster-larvae-6221663"

REPO_ROOT = Path(__file__).resolve().parents[3]
STUDY_DIR = REPO_ROOT / "data" / "studies" / STUDY_ID
SAMPLESHEET = STUDY_DIR / "samplesheet.tsv"
PROVENANCE = STUDY_DIR / "pride_provenance.json"

# (MSstats Condition, BioReplicate) -> hatchery proteomic sample number.
SKYLINE_REPLICATES = {("23C", "1"): 1, ("29C", "1"): 2, ("23C", "2"): 13, ("29C", "2"): 14}
# Exposure day 0 is 6 Sep 2015, when the 29 C silos were set (hatchery Notes sheet).
DAY0 = date(2015, 9, 6)
# Each pool was injected four times, one precursor m/z range per injection
# (PRIDE protocol: 400-525, 535-650, 650-775, 775-900). Which file is which range
# is not recorded, but every pool has one file at each of four MSstats-report row
# counts (the Skyline document is shared), so the row count groups the four
# injections across pools. injection_group numbers those groups by row count;
# it is not the m/z order.
CONDITION = {"23C": "control", "29C": "treatment"}


def fetch(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "AREE-curation"})
    with urllib.request.urlopen(req, timeout=600) as resp:
        return resp.read()


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def hatchery_key(xlsx: bytes) -> dict[int, tuple[date, str]]:
    ws = openpyxl.load_workbook(io.BytesIO(xlsx), data_only=True)["Sheet1"]
    key = {}
    for row in ws.iter_rows(values_only=True):
        cells = [c for c in row if c is not None]
        if len(cells) >= 3 and isinstance(cells[1], datetime) and str(cells[0]).isdigit():
            key[int(cells[0])] = (cells[1].date(), str(cells[2]))
    return key


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--download", type=Path, required=True, help="staging directory for fetched files")
    args = ap.parse_args()
    stage = args.download
    stage.mkdir(parents=True, exist_ok=True)

    project_raw = fetch(PRIDE_API)
    files_raw = fetch(f"{PRIDE_API}/files/all")
    record_raw = fetch(ZENODO_API)
    for name, data in (("pride_project.json", project_raw), ("pride_files.json", files_raw),
                       ("zenodo_record.json", record_raw)):
        (stage / name).write_bytes(data)
    project, files, record = json.loads(project_raw), json.loads(files_raw), json.loads(record_raw)

    (zfile,) = record["files"]
    md5 = zfile["checksum"].removeprefix("md5:")
    zip_path = stage / Path(zfile["key"]).name
    if not zip_path.exists() or hashlib.md5(zip_path.read_bytes()).hexdigest() != md5:
        zip_path.write_bytes(fetch(zfile["links"]["self"]))
    zip_bytes = zip_path.read_bytes()
    if hashlib.md5(zip_bytes).hexdigest() != md5:
        sys.exit(f"Zenodo zip md5 does not match the Zenodo API ({md5})")
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
        report_zip = zf.read(f"{ZIP_ROOT}/data/20190403-MSstats-report.csv.zip")
        protcomp_raw = zf.read(f"{ZIP_ROOT}/data/20190403-2015Cgseed-protcomp.csv")
        xlsx_raw = zf.read(f"{ZIP_ROOT}/data/2015X15P8LargescaleSeed29CtoAmbient.xlsx")
        msstats_script = zf.read(f"{ZIP_ROOT}/scripts/01-MSstats-Cgseed-DIA.R").decode()
    with zipfile.ZipFile(io.BytesIO(report_zip)) as rz:
        report_raw = rz.read("20190403-MSstats-report.csv")

    # Run file -> (Condition, BioReplicate), and rows per file (identifies the m/z range).
    runs: dict[str, tuple[str, str]] = {}
    rows_per_file: Counter = Counter()
    for r in csv.DictReader(io.StringIO(report_raw.decode())):
        f = r["File Name"].removesuffix(".mzML")
        cb = (r["Condition"], r["BioReplicate"])
        if runs.setdefault(f, cb) != cb:
            sys.exit(f"{f} has more than one Condition/BioReplicate in the MSstats report")
        rows_per_file[f] += 1
    by_name = {f["fileName"]: f for f in files}
    raw_names = {n.removesuffix(".raw") for n in by_name if n.endswith(".raw")}
    if set(runs) != raw_names:
        sys.exit(f"MSstats report and PRIDE .raw files differ: {sorted(set(runs) ^ raw_names)}")
    range_sizes = sorted(set(rows_per_file.values()))
    if len(range_sizes) != 4:
        sys.exit(f"expected four injection ranges per pool, got row counts {range_sizes}")

    key = hatchery_key(xlsx_raw)
    pools = {}
    for (cond, biorep), number in SKYLINE_REPLICATES.items():
        sampled, contents = key[number]
        silo = re.search(r"Silo (\d+)", contents).group(1)
        expected = {"23C": ("23C control seed", "9"), "29C": ("29C control seed", "8")}[cond]
        if not contents.startswith(expected[0]) or silo != expected[1]:
            sys.exit(f"hatchery sample {number} is {contents!r}, not the {cond} control silo")
        pools[(cond, biorep)] = {"hatchery_sample": number, "silo": silo, "date": sampled,
                                 "day": (sampled - DAY0).days}
    if {p["day"] for p in pools.values()} != {5, 8}:
        sys.exit(f"expected sampling on exposure days 5 and 8, got {sorted(p['day'] for p in pools.values())}")

    out = []
    for run, (cond, biorep) in runs.items():
        pool = pools[(cond, biorep)]
        raw_file = f"{run}.raw"
        out.append({
            "run_id": run,
            "pool_id": f"Oysterseed_{pool['hatchery_sample']:02d}",
            "condition": CONDITION[cond],
            "temperature_C": cond.removesuffix("C"),
            "silo": pool["silo"],
            "hatchery_sample": pool["hatchery_sample"],
            "date_sampled": pool["date"].isoformat(),
            "day_of_exposure": pool["day"],
            "msstats_bioreplicate": biorep,
            "injection_group": range_sizes.index(rows_per_file[run]) + 1,
            "sample_unit": "pool of dry settled seed (~500 uL) from one silo",
            "raw_file": raw_file,
            "raw_sha1": by_name[raw_file]["checksum"],
            "raw_bytes": by_name[raw_file]["fileSizeBytes"],
            "elib_file": f"{run}.mzML.elib",
        })
    out.sort(key=lambda r: (r["day_of_exposure"], r["temperature_C"], r["injection_group"]))

    pool_counts = Counter(r["pool_id"] for r in out)
    if len(pool_counts) != 4 or set(pool_counts.values()) != {4}:
        sys.exit(f"expected 4 pools x 4 injections, got {dict(pool_counts)}")

    with open(SAMPLESHEET, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    comp = list(csv.DictReader(io.StringIO(protcomp_raw.decode())))
    adj = [float(r["adj.pvalue"]) for r in comp if r["adj.pvalue"] not in ("", "NA")]
    contrast = re.search(r"matrix\(c\((-?1), *(-?1)\)", msstats_script)
    provenance = {
        "study_id": STUDY_ID,
        "accession": ACCESSION,
        "repository": "PRIDE",
        "license": project["license"],
        "submission_type": project["submissionType"],
        "date_retrieved": date.today().isoformat(),
        "generated_by": f"data/studies/{STUDY_ID}/build_samplesheet.py",
        "source_urls": [PRIDE_API, f"{PRIDE_API}/files/all", ZENODO_API],
        "pride_reference": project["references"][0]["referenceLine"],
        "pride_instrument": [i["name"] for i in project["instruments"]],
        "pride_files_by_category": dict(sorted(Counter(f["fileCategory"]["value"] for f in files).items())),
        "input_sha256": {
            Path(zfile["key"]).name: sha256(zip_bytes),
            "zip:data/20190403-MSstats-report.csv": sha256(report_raw),
            "zip:data/20190403-2015Cgseed-protcomp.csv": sha256(protcomp_raw),
            "zip:data/2015X15P8LargescaleSeed29CtoAmbient.xlsx": sha256(xlsx_raw),
        },
        "zenodo": {
            "doi": record["doi"],
            "version": record["metadata"].get("version"),
            "license": record["metadata"]["license"]["id"],
            "file": zfile["key"],
            "md5": md5,
            "size_bytes": zfile["size"],
        },
        "n_runs": len(out),
        "n_pools": len(pool_counts),
        "injections_per_pool": 4,
        "pools": {f"Oysterseed_{p['hatchery_sample']:02d}": {
            "condition": c, "silo": p["silo"], "date": p["date"].isoformat(), "day_of_exposure": p["day"]}
            for (c, _), p in sorted(pools.items(), key=lambda kv: kv[1]["hatchery_sample"])},
        "silos_by_temperature_C": {"23": ["9"], "29": ["8"]},
        "msstats_result": {
            "n_proteins": len(comp),
            "label": sorted({r["Label"] for r in comp}),
            "contrast_vector_23C_29C": [int(contrast.group(1)), int(contrast.group(2))],
            "n_adj_p_below_0.05": sum(a < 0.05 for a in adj),
            "df_counts": dict(Counter(r["DF"] for r in comp).most_common()),
        },
        "outputs": {f"data/studies/{STUDY_ID}/samplesheet.tsv": sha256(SAMPLESHEET.read_bytes())},
    }
    PROVENANCE.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"wrote {SAMPLESHEET.relative_to(REPO_ROOT)} ({len(out)} runs) and {PROVENANCE.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
