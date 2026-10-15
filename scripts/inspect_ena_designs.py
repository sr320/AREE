#!/usr/bin/env python3
"""Recover the deposited design of candidate BioProjects from ENA metadata.

The second step of the screen in docs/candidate_studies_mussels_clams.md. For
each project: the species actually recorded on each run (catches mislabelled
Mytilus deposits), layout and depth, the study description and any PubMed
links, and every sample attribute that varies, with group sizes. Downloads no
sequence data.

Usage (repo root, aree installed): python scripts/inspect_ena_designs.py <out.txt> ACC [ACC ...]
"""
from __future__ import annotations

import statistics
import sys
import time
import xml.etree.ElementTree as ET
from collections import Counter

from intake.ena_samplesheet import ENA_XML, _get, fetch_run_report, fetch_sample_attributes

# Attributes that vary between samples but never describe the design.
SKIP = {
    "_alias", "ENA-FIRST-PUBLIC", "ENA-LAST-UPDATE", "INSDC first public", "INSDC last update",
    "INSDC center name", "INSDC status", "SRA accession", "INSDC secondary accession",
    "Submitter Id", "ENA-CHECKLIST", "collection_date", "lat_lon", "geo_loc_name",
    "BioSampleModel", "organism", "isolate_name_alias",
}


def study_info(accession: str) -> tuple[str, list[str], str]:
    try:
        root = ET.fromstring(_get(f"{ENA_XML}/{accession}"))  # noqa: S314 - EBI XML
    except (OSError, ET.ParseError) as exc:
        return "", [], f"(study XML failed: {exc})"
    title = root.findtext(".//TITLE") or ""
    description = (root.findtext(".//DESCRIPTION") or "").strip().replace("\n", " ")
    pubmed = [x.findtext("ID") for x in root.iter("XREF_LINK") if (x.findtext("DB") or "").upper() == "PUBMED"]
    return title, pubmed, description


def describe(accession: str) -> list[str]:
    try:
        runs = fetch_run_report(accession, ["scientific_name", "tax_id", "sample_title",
                                            "experiment_title", "library_strategy"])
    except OSError as exc:
        return [f"\n### {accession}: run report failed: {exc}"]
    rnaseq = [r for r in runs if r.get("library_strategy") == "RNA-Seq"]
    samples = sorted({r["sample_accession"] for r in rnaseq})
    attrs = fetch_sample_attributes(samples) if samples else {}
    title, pubmed, description = study_info(accession)
    species = Counter(f"{r['scientific_name']}({r['tax_id']})" for r in rnaseq)
    layouts = Counter(r.get("library_layout") for r in rnaseq)
    reads = [int(r["read_count"]) for r in rnaseq if r.get("read_count", "").isdigit()]
    median_m = statistics.median(reads) / 1e6 if reads else 0

    lines = [
        f"\n### {accession} | {title}",
        f"runs={len(runs)} rnaseq_runs={len(rnaseq)} samples={len(samples)} layout={dict(layouts)} "
        f"median_reads={median_m:.1f}M species={dict(species)} pubmed={pubmed}",
        f"desc: {description[:600]}",
    ]
    keys = Counter(k for s in samples for k in (attrs.get(s) or {}) if k not in SKIP)
    for key, _ in keys.most_common():
        values = Counter(str((attrs.get(s) or {}).get(key)) for s in samples)
        if 1 < len(values) <= 40:
            lines.append(f"  {key}: {dict(values.most_common(14))}")
    lines.append(f"  aliases: {[(attrs.get(s) or {}).get('_alias') for s in samples][:8]}")
    lines.append(f"  sample_titles: {sorted({r.get('sample_title', '') for r in rnaseq})[:8]}")
    return lines


def main(out_path: str, accessions: list[str]) -> None:
    with open(out_path, "w") as out:
        for accession in accessions:
            out.write("\n".join(describe(accession)) + "\n")
            out.flush()
            time.sleep(1)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
