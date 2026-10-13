#!/usr/bin/env python3
"""List every public RNA-Seq read run in ENA for the mussel and clam taxa, by project.

The first step of the screen in docs/candidate_studies_mussels_clams.md. Reads
only ENA's portal API; downloads no sequence data.

Usage: python scripts/list_ena_rnaseq_projects.py <out_dir>
Writes runs_<slug>.tsv (one row per run) and projects.tsv (one row per project).
"""
from __future__ import annotations

import csv
import io
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

TAXA = {"mgall": 29158, "medul": 6550, "mtros": 6551, "mcali": 6549, "mmerc": 6596, "rphil": 129788}
FIELDS = (
    "run_accession,study_accession,secondary_study_accession,sample_accession,sample_alias,"
    "sample_title,experiment_title,run_alias,library_layout,library_source,library_selection,"
    "instrument_model,read_count,base_count,fastq_bytes,first_public,study_title,scientific_name,tax_id"
)
API = "https://www.ebi.ac.uk/ena/portal/api/search?"


def get(url: str, tries: int = 5) -> str:
    """ENA returns intermittent 500s for large tax_tree queries; retry with backoff."""
    for attempt in range(tries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "aree-curation"})
            with urllib.request.urlopen(request, timeout=180) as response:  # noqa: S310 - fixed EBI host
                return response.read().decode()
        except OSError as exc:
            error = exc
            time.sleep(5 * (attempt + 1))
    raise error


def fastq_gb(run: dict) -> float:
    return sum(int(x) for x in run["fastq_bytes"].split(";") if x) / 1e9


def main(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    projects: dict[str, dict] = {}
    for slug, taxid in TAXA.items():
        query = f'tax_tree({taxid}) AND library_strategy="RNA-Seq"'
        text = get(API + urllib.parse.urlencode(
            {"result": "read_run", "query": query, "fields": FIELDS, "format": "tsv", "limit": "0"}))
        (out / f"runs_{slug}.tsv").write_text(text)
        rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
        by_study = defaultdict(list)
        for run in rows:
            by_study[run["study_accession"]].append(run)
        for study, runs in by_study.items():
            p = projects.setdefault(study, {
                "taxa": set(), "runs": 0, "samples": set(), "gb": 0.0, "layouts": set(),
                "instruments": set(), "first_public": "9999",
                "secondary": runs[0]["secondary_study_accession"], "title": runs[0]["study_title"],
            })
            p["taxa"].add(slug)
            p["runs"] += len(runs)
            p["samples"].update(r["sample_accession"] for r in runs)
            p["gb"] += sum(fastq_gb(r) for r in runs)
            p["layouts"].update(r["library_layout"] for r in runs)
            p["instruments"].update(r["instrument_model"] for r in runs)
            p["first_public"] = min(p["first_public"], *(r["first_public"] for r in runs))
        print(slug, taxid, len(rows), "runs", len(by_study), "projects", file=sys.stderr)
        time.sleep(1)

    with open(out / "projects.tsv", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["study_accession", "secondary", "taxa", "runs", "samples", "fastq_gb",
                         "layout", "instrument", "first_public", "title"])
        for study, p in sorted(projects.items(), key=lambda kv: (sorted(kv[1]["taxa"]), -kv[1]["runs"])):
            writer.writerow([study, p["secondary"], ",".join(sorted(p["taxa"])), p["runs"],
                             len(p["samples"]), f"{p['gb']:.1f}", ",".join(sorted(p["layouts"])),
                             ",".join(sorted(p["instruments"]))[:60], p["first_public"], p["title"]])
    print(len(projects), "projects", file=sys.stderr)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
