"""Build a reproducible sample sheet and FASTQ manifest from an ENA/SRA BioProject.

A raw-reanalysis study starts with a question the deposited metadata can answer
but a paper abstract usually cannot: which runs belong to which experimental
group, and how many replicates each group actually has. Hand-transcribing that
from a supplementary methods table is where curation errors enter, and it is
also how an unreplicated design slips through unnoticed (see
docs/candidate_studies.md on PRJNA623063).

This module reads the authoritative record instead — the ENA read-run report
plus each sample's attribute block — and emits:

* ``samplesheet.tsv``      — one row per run, with the sample attributes that
                             define the design, in the column contract
                             ``workflows/rnaseq`` expects;
* ``fastq_manifest.tsv``   — FASTQ URLs with ENA's own MD5s and byte counts, so
                             a download can be verified rather than trusted;
* ``ena_provenance.json``  — the queries used, the date, and a checksum of each
                             emitted file.

Nothing here downloads sequence data. It produces the plan for doing so.
"""
from __future__ import annotations

import json
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path

from common import REPO_ROOT, sha256_file

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_XML = "https://www.ebi.ac.uk/ena/browser/api/xml"

RUN_FIELDS = [
    "run_accession", "sample_accession", "sample_alias", "experiment_accession",
    "library_strategy", "library_layout", "library_selection", "instrument_model",
    "read_count", "base_count", "fastq_ftp", "fastq_md5", "fastq_bytes",
]

TIMEOUT = 120


class ENAError(RuntimeError):
    """A query returned nothing usable, or the project is not what was expected."""


def _get(url: str) -> str:
    with urllib.request.urlopen(url, timeout=TIMEOUT) as fh:  # noqa: S310 - fixed EBI host
        return fh.read().decode("utf-8")


def _run_fields(extra: list[str] | None = None) -> list[str]:
    return list(dict.fromkeys(RUN_FIELDS + (extra or [])))


def fetch_run_report(bioproject: str, extra_fields: list[str] | None = None) -> list[dict]:
    """One record per sequencing run, from ENA's read_run report."""
    query = urllib.parse.urlencode(
        {"accession": bioproject, "result": "read_run",
         "fields": ",".join(_run_fields(extra_fields)), "format": "tsv"}
    )
    text = _get(f"{ENA_PORTAL}?{query}")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise ENAError(
            f"ENA returned no runs for {bioproject}. Either the accession is wrong, "
            "or the project has no public reads (some BioProjects register metadata "
            "only)."
        )
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]


def fetch_sample_attributes(sample_accessions: list[str]) -> dict:
    """Attribute block per sample: the experimental factors, as deposited."""
    out: dict[str, dict] = {}
    # ENA accepts comma-separated accessions but the URL has a practical limit.
    for i in range(0, len(sample_accessions), 40):
        batch = sample_accessions[i : i + 40]
        root = ET.fromstring(_get(f"{ENA_XML}/{','.join(batch)}"))
        for sample in root.iter("SAMPLE"):
            primary = sample.findtext(".//PRIMARY_ID")
            ext = [e.text for e in sample.iter("EXTERNAL_ID") if e.text]
            attrs = {
                a.findtext("TAG"): a.findtext("VALUE")
                for a in sample.iter("SAMPLE_ATTRIBUTE")
                if a.findtext("TAG")
            }
            attrs["_alias"] = sample.get("alias")
            for key in filter(None, [primary, *ext]):
                out[key] = attrs
    return out


def _slug(value) -> str:
    return "".join(c if c.isalnum() else "_" for c in str(value or "NA")).strip("_")


def parse_run_labels(runs: list[dict], label_field: str, label_pattern: str) -> dict[str, dict]:
    """Per-run design attributes parsed from a run-level field.

    For projects that deposit every library under one BioSample, the design
    lives only in a run- or experiment-level string (e.g. experiment_title
    "... sequencing: RR52-R1"). The regex must match every run and its named
    groups become attributes; the matched text becomes the run's sample_id.
    A run that does not match is an error, never a silent default.
    """
    rx = re.compile(label_pattern)
    if not rx.groupindex:
        raise ENAError("--label-pattern must define at least one named group, e.g. (?P<family>\\d+).")
    out, unmatched = {}, []
    for r in runs:
        m = rx.search(r.get(label_field) or "")
        if not m:
            unmatched.append(f"{r['run_accession']} ({label_field}={r.get(label_field)!r})")
            continue
        out[r["run_accession"]] = {**m.groupdict(), "_label": m.group(0)}
    if unmatched:
        raise ENAError(f"{len(unmatched)} run(s) do not match --label-pattern: {unmatched[:5]}")
    labels = [v["_label"] for v in out.values()]
    if len(set(labels)) != len(labels):
        raise ENAError("--label-pattern matched the same label on more than one run; "
                       "labels must identify libraries uniquely.")
    return out


def write_methylkit_samplesheet(rows: list[dict], path: Path, control_level: str, treatment_level: str) -> dict:
    """methylKit sample sheet (sample_id,treatment) keyed by run accession.

    workflows/methylation pairs FASTQ by file prefix (the run accession), so
    sample_id here is the run accession, not the curated label. Runs in any
    other condition are left out rather than coded.
    """
    coded = {control_level: 0, treatment_level: 1}
    missing = [lvl for lvl in coded if lvl not in {r["condition"] for r in rows}]
    if missing:
        raise ENAError(f"condition level(s) {missing} not found in the sample sheet.")
    kept = sorted((r for r in rows if r["condition"] in coded), key=lambda r: r["run_accession"])
    path.write_text("sample_id,treatment\n" + "".join(
        f"{r['run_accession']},{coded[r['condition']]}\n" for r in kept))
    return {"control_level": control_level, "treatment_level": treatment_level,
            "n_control": sum(r["condition"] == control_level for r in kept),
            "n_treatment": sum(r["condition"] == treatment_level for r in kept)}


def build_samplesheet(
    bioproject: str,
    *,
    study_id: str,
    out_dir: Path,
    condition_attributes: list[str],
    extra_attributes: list[str] | None = None,
    label_field: str | None = None,
    label_pattern: str | None = None,
    methylkit_levels: tuple[str, str] | None = None,
) -> dict:
    """Write samplesheet + FASTQ manifest + provenance for one BioProject.

    `condition_attributes` names the sample attributes that jointly define an
    experimental group; their slugified concatenation becomes the `condition`
    column the DESeq2 step models. Naming them explicitly — rather than
    inferring groups from sample aliases — keeps the design decision visible in
    the curation record instead of buried in a string-splitting heuristic.

    When `label_field` and `label_pattern` are given, the named groups of the
    pattern (see `parse_run_labels`) are also available as attributes; they
    take precedence over sample attributes of the same name. `methylkit_levels`
    = (control, treatment) additionally writes methylkit_samplesheet.csv.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    extra_attributes = extra_attributes or []
    if bool(label_field) != bool(label_pattern):
        raise ENAError("--label-field and --label-pattern must be given together.")

    runs = fetch_run_report(bioproject, [label_field] if label_field else None)
    attrs = fetch_sample_attributes(sorted({r["sample_accession"] for r in runs}))
    run_labels = parse_run_labels(runs, label_field, label_pattern) if label_field else {}

    strategies = {r["library_strategy"] for r in runs}
    layouts = {r["library_layout"] for r in runs}

    all_attr_cols = list(dict.fromkeys(condition_attributes + extra_attributes))
    sheet_cols = ["sample_id", "run_accession", "sample_accession", "sample_alias",
                  "condition", *[_slug(c).lower() for c in all_attr_cols],
                  "library_layout", "instrument_model", "read_count"]

    rows, manifest_rows, missing = [], [], []
    for r in sorted(runs, key=lambda x: x["run_accession"]):
        a = {**attrs.get(r["sample_accession"], {}), **run_labels.get(r["run_accession"], {})}
        values = [a.get(c) for c in condition_attributes]
        if any(v is None for v in values):
            missing.append(f"{r['run_accession']} ({r['sample_accession']})")
        condition = "_".join(_slug(v) for v in values)

        rows.append([
            a.get("_label") or r["sample_alias"] or r["run_accession"], r["run_accession"],
            r["sample_accession"], r["sample_alias"], condition,
            *[str(a.get(c) or "") for c in all_attr_cols],
            r["library_layout"], r["instrument_model"], r["read_count"],
        ])

        urls = [u for u in (r.get("fastq_ftp") or "").split(";") if u]
        md5s = [m for m in (r.get("fastq_md5") or "").split(";") if m]
        sizes = [b for b in (r.get("fastq_bytes") or "").split(";") if b]
        for idx, url in enumerate(urls):
            manifest_rows.append([
                r["run_accession"], str(idx + 1), f"ftp://{url}",
                md5s[idx] if idx < len(md5s) else "",
                sizes[idx] if idx < len(sizes) else "",
            ])

    if missing:
        raise ENAError(
            f"{len(missing)} run(s) lack one of the condition attributes "
            f"{condition_attributes}: {missing[:5]}... "
            "Check the attribute names against the deposited sample records — "
            "they are submitter-defined and vary between projects."
        )

    sheet_path = out_dir / "samplesheet.tsv"
    sheet_path.write_text(
        "\t".join(sheet_cols) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n"
    )
    manifest_path = out_dir / "fastq_manifest.tsv"
    manifest_path.write_text(
        "run_accession\tmate\turl\tmd5\tbytes\n"
        + "\n".join("\t".join(r) for r in manifest_rows) + "\n"
    )

    group_counts: dict[str, int] = {}
    for r in rows:
        group_counts[r[4]] = group_counts.get(r[4], 0) + 1

    total_bytes = sum(int(m[4]) for m in manifest_rows if m[4])

    output_paths = [sheet_path, manifest_path]
    methylkit = None
    if methylkit_levels:
        mk_path = out_dir / "methylkit_samplesheet.csv"
        methylkit = write_methylkit_samplesheet(
            [{"run_accession": r[1], "condition": r[4]} for r in rows], mk_path, *methylkit_levels)
        output_paths.append(mk_path)
    provenance = {
        "study_id": study_id,
        "bioproject": bioproject,
        "date_generated": date.today().isoformat(),
        "generated_by": "src/intake/ena_samplesheet.py",
        "queries": {
            "run_report": f"{ENA_PORTAL}?accession={bioproject}&result=read_run&fields={','.join(_run_fields([label_field] if label_field else None))}",
            "sample_xml": f"{ENA_XML}/<sample_accessions>",
        },
        "condition_attributes": condition_attributes,
        "n_runs": len(rows),
        "n_samples": len({r[2] for r in rows}),
        "library_strategies": sorted(strategies),
        "library_layouts": sorted(layouts),
        "group_sizes": dict(sorted(group_counts.items())),
        "min_group_size": min(group_counts.values()) if group_counts else 0,
        "total_fastq_bytes": total_bytes,
        "outputs": {
            str(p.relative_to(REPO_ROOT)) if p.is_relative_to(REPO_ROOT) else str(p): sha256_file(p)
            for p in output_paths
        },
    }
    if label_field:
        # Recorded only when used, so earlier provenance files stay reproducible.
        provenance["run_labels"] = {"field": label_field, "pattern": label_pattern,
                                    "labels": {k: v["_label"] for k, v in sorted(run_labels.items())}}
    if methylkit:
        provenance["methylkit_samplesheet"] = methylkit
    prov_path = out_dir / "ena_provenance.json"
    prov_path.write_text(json.dumps(provenance, indent=2) + "\n")
    provenance["provenance_file"] = str(prov_path)
    return provenance
