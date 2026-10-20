"""Build the sample sheet and provenance for WANAMAKER2020_HEAT_SEED (PXD013262).

PRIDE PXD013262 (CC0) deposits 26 LC-MS/MS runs (.raw, .mzXML, Comet
.pep.xml) and a short metadata CSV giving day and temperature per run. It does
not say which rearing silo a run came from, or that runs come in technical
pairs. The authors' supplementary package on figshare
(10.6084/m9.figshare.12436598, CC BY 4.0) does: its Data/Sample_metadata.csv
names every run's pool, silo and technical replicate, and Data/Abacus_output.tsv
is the protein x run abundance table (adjusted NSAF) the paper analysed.

    python -I data/studies/WANAMAKER2020_HEAT_SEED/build_samplesheet.py --download STAGING_DIR

fetches the PRIDE project and file records, the PRIDE metadata CSV and the
figshare zip (about 40 MB, md5-checked against the figshare API) into
STAGING_DIR, then writes

  data/studies/WANAMAKER2020_HEAT_SEED/samplesheet.tsv      one row per MS run
  data/studies/WANAMAKER2020_HEAT_SEED/pride_provenance.json

Nothing else is committed. The script stops if the two metadata sources
disagree on any run's day or temperature, or if the design is not one silo per
temperature with two technical runs per pool, which is the finding the study's
QC decision rests on.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
import urllib.request
import zipfile
from collections import Counter
from datetime import date
from pathlib import Path

ACCESSION = "PXD013262"
STUDY_ID = "WANAMAKER2020_HEAT_SEED"
PRIDE_API = f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{ACCESSION}"
PRIDE_FTP = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/06/PXD013262"
PRIDE_METADATA = "oyster_tempXtime_proteomics_metadata.csv"
FIGSHARE_ARTICLE = 12436598
FIGSHARE_API = f"https://api.figshare.com/v2/articles/{FIGSHARE_ARTICLE}"
ZIP_ROOT = "shellytrigg-paper-OysterSeed-TimeXTemp-73d4ad3"

REPO_ROOT = Path(__file__).resolve().parents[3]
STUDY_DIR = REPO_ROOT / "data" / "studies" / STUDY_ID
SAMPLESHEET = STUDY_DIR / "samplesheet.tsv"
PROVENANCE = STUDY_DIR / "pride_provenance.json"

SILO_TEMPERATURE = {"3": 23, "9": 29}
CONDITION = {16: "baseline", 23: "control", 29: "treatment"}


def fetch(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "AREE-curation"})
    with urllib.request.urlopen(req, timeout=300) as resp:
        return resp.read()


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_csv(data: bytes) -> list[dict]:
    # Both metadata files were saved with classic-Mac CR line endings.
    text = data.decode("utf-8-sig").replace("\r\n", "\n").replace("\r", "\n")
    return list(csv.DictReader(io.StringIO(text)))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--download", type=Path, required=True, help="staging directory for fetched files")
    args = ap.parse_args()
    stage = args.download
    stage.mkdir(parents=True, exist_ok=True)

    project_raw = fetch(PRIDE_API)
    files_raw = fetch(f"{PRIDE_API}/files/all")
    meta_raw = fetch(f"{PRIDE_FTP}/{PRIDE_METADATA}")
    article_raw = fetch(FIGSHARE_API)
    for name, data in (("pride_project.json", project_raw), ("pride_files.json", files_raw),
                       (PRIDE_METADATA, meta_raw), ("figshare_article.json", article_raw)):
        (stage / name).write_bytes(data)
    project, files, article = json.loads(project_raw), json.loads(files_raw), json.loads(article_raw)

    (fs_file,) = article["files"]
    zip_path = stage / fs_file["name"]
    if not zip_path.exists() or hashlib.md5(zip_path.read_bytes()).hexdigest() != fs_file["computed_md5"]:
        zip_path.write_bytes(fetch(fs_file["download_url"]))
    zip_bytes = zip_path.read_bytes()
    if hashlib.md5(zip_bytes).hexdigest() != fs_file["computed_md5"]:
        sys.exit(f"figshare zip md5 does not match the figshare API ({fs_file['computed_md5']})")
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
        fs_meta_raw = zf.read(f"{ZIP_ROOT}/Data/Sample_metadata.csv")
        abacus_raw = zf.read(f"{ZIP_ROOT}/Data/Abacus_output.tsv")
        members = set(zf.namelist())

    by_name = {f["fileName"]: f for f in files}
    categories = Counter(f["fileCategory"]["value"] for f in files)

    # PRIDE metadata: one row per deposited file; keep the .raw rows.
    pride_runs = {r["SampleID"]: r for r in read_csv(meta_raw) if r["filename"].endswith(".raw")}
    fs_runs = {r["SampleID"]: r for r in read_csv(fs_meta_raw)}
    if set(pride_runs) != set(fs_runs):
        sys.exit(f"PRIDE and figshare list different runs: {sorted(set(pride_runs) ^ set(fs_runs))}")

    abacus_header = abacus_raw.decode().split("\n", 1)[0].split("\t")
    abacus_runs = {c.removeprefix("20161205_SAMPLE_").removesuffix("_ADJNSAF")
                   for c in abacus_header if c.endswith("_ADJNSAF")}
    if abacus_runs != {r.upper() for r in fs_runs}:
        sys.exit("Abacus_output.tsv adjusted-NSAF columns do not match the 26 deposited runs")
    n_proteins = abacus_raw.decode().rstrip("\n").count("\n")

    rows = []
    for run_id, fs in fs_runs.items():
        pr = pride_runs[run_id]
        contents = fs["Contents"]                      # e.g. Silo3_23C_21dpf_Rep1, 19dpf_Rep2
        rep = int(contents.rsplit("_Rep", 1)[1])
        if contents.startswith("Silo"):
            silo, temp_label, dpf_label = contents.split("_")[:3]
            silo = silo.removeprefix("Silo")
            temp = int(temp_label.removesuffix("C"))
            if SILO_TEMPERATURE[silo] != temp:
                sys.exit(f"run {run_id}: silo {silo} is not the {temp} C silo")
        else:
            silo, temp, dpf_label = "", 16, contents.split("_")[0]
        dpf = int(dpf_label.removesuffix("dpf"))
        if int(pr["temp"]) != temp:
            sys.exit(f"run {run_id}: PRIDE says {pr['temp']} C, figshare says {temp} C")
        # Day 0 is the 19 dpf baseline; exposure day d is d + 18 dpf (paper, Methods).
        if dpf != (19 if int(pr["day"]) == 0 else int(pr["day"]) + 18):
            sys.exit(f"run {run_id}: PRIDE day {pr['day']} does not match {dpf} dpf")
        raw = f"20161205_{'sample' if run_id == '1' else 'Sample'}_{run_id}"
        raw_file = pr["filename"]
        stem = raw_file.removesuffix(".raw")
        peak = next(n for n in (f"{stem}.mzXML", f"{raw}.mzXML") if n in by_name)
        search = next(n for n in (f"{stem}.pep.xml", f"{raw}.pep.xml") if n in by_name)
        rows.append({
            "run_id": run_id,
            "pool_id": fs["SampleName"],
            "condition": CONDITION[temp],
            "temperature_C": temp,
            "silo": silo,
            "day_of_exposure": int(pr["day"]),
            "dpf": dpf,
            "technical_replicate": rep,
            "uwpr_alternate_id": fs["AlternateID"],
            "sample_unit": ("pool of ~12,500 larvae/seed from one silo" if silo
                            else "pool of ~12,500 pediveliger larvae before transfer to silos"),
            "raw_file": raw_file,
            "raw_sha1": by_name[raw_file]["checksum"],
            "raw_bytes": by_name[raw_file]["fileSizeBytes"],
            "peak_file": peak,
            "search_file": search,
        })
    rows.sort(key=lambda r: (r["day_of_exposure"], r["temperature_C"], r["technical_replicate"]))

    # The design facts the QC decision rests on.
    pools = Counter(r["pool_id"] for r in rows)
    if set(pools.values()) != {2} or len(pools) != 13:
        sys.exit(f"expected 13 pools x 2 technical runs, got {dict(pools)}")
    silos_by_temp = {t: sorted({r["silo"] for r in rows if r["temperature_C"] == t}) for t in (23, 29)}
    if silos_by_temp != {23: ["3"], 29: ["9"]}:
        sys.exit(f"expected one silo per temperature, got {silos_by_temp}")
    pools_by_temp_day = Counter((r["temperature_C"], r["day_of_exposure"]) for r in rows if r["silo"])
    if set(pools_by_temp_day.values()) != {2} or len(pools_by_temp_day) != 12:
        sys.exit("expected one pool (two runs) per temperature x day")

    with open(SAMPLESHEET, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    provenance = {
        "study_id": STUDY_ID,
        "accession": ACCESSION,
        "repository": "PRIDE",
        "license": project["license"],
        "submission_type": project["submissionType"],
        "date_retrieved": date.today().isoformat(),
        "generated_by": f"data/studies/{STUDY_ID}/build_samplesheet.py",
        "source_urls": [PRIDE_API, f"{PRIDE_API}/files/all", f"{PRIDE_FTP}/{PRIDE_METADATA}", FIGSHARE_API],
        "pride_reference": project["references"][0]["referenceLine"],
        "pride_instrument": [i["name"] for i in project["instruments"]],
        "pride_files_by_category": dict(sorted(categories.items())),
        "search_database_deposited_in_pride": "Cg-Giga_cont_AA.fa" in by_name,
        "search_database_in_figshare": f"{ZIP_ROOT}/Data/Cg-Giga_cont_AA.fa" in members,
        "input_sha256": {
            PRIDE_METADATA: sha256(meta_raw),
            fs_file["name"]: sha256(zip_bytes),
            f"{fs_file['name']}:Data/Sample_metadata.csv": sha256(fs_meta_raw),
            f"{fs_file['name']}:Data/Abacus_output.tsv": sha256(abacus_raw),
        },
        "figshare": {
            "doi": article["doi"],
            "version": article["version"],
            "license": article["license"]["name"],
            "file": fs_file["name"],
            "md5": fs_file["computed_md5"],
            "size_bytes": fs_file["size"],
        },
        "n_runs": len(rows),
        "n_pools": len(pools),
        "technical_runs_per_pool": 2,
        "silos_by_temperature_C": {str(t): s for t, s in silos_by_temp.items()},
        "pools_per_temperature_per_day": 1,
        "days_of_exposure": sorted({r["day_of_exposure"] for r in rows if r["silo"]}),
        "abacus_output": {"n_proteins": n_proteins, "n_adjnsaf_columns": len(abacus_runs)},
        "outputs": {f"data/studies/{STUDY_ID}/samplesheet.tsv": sha256(SAMPLESHEET.read_bytes())},
    }
    PROVENANCE.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"wrote {SAMPLESHEET.relative_to(REPO_ROOT)} ({len(rows)} runs) and {PROVENANCE.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
