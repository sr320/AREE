"""Convert a published supplementary differential-expression table into AREE result files.

Real studies rarely publish a DESeq2 output in the exact shape the harmonizers
expect. This module performs the *minimum* mechanical reshaping needed —
selecting and renaming columns, dropping trailing non-data rows — and records
what it did, so that the transformation between the published artifact and the
harmonized evidence is inspectable rather than a manual copy-paste.

It deliberately does NOT:

* alter identifiers (decoration such as ``gene-LOC123|LOC123`` is preserved
  verbatim; `harmonize.identifiers` handles lookup),
* impute missing statistics (a published table that reports only an adjusted
  p-value yields a results file with no raw p-value and no standard error),
* apply any additional filtering beyond what the source authors applied.
"""
from __future__ import annotations

import json
from datetime import date
from pathlib import Path

import pandas as pd

from common import REPO_ROOT, sha256_file

# AREE's RNA-seq harmonizer contract (see src/harmonize/rnaseq.py).
RNASEQ_COLUMNS = ["gene_id", "log2FoldChange", "lfcSE", "pvalue", "padj"]


def convert_de_table(
    frame: pd.DataFrame,
    column_map: dict,
    *,
    out_path: Path,
) -> dict:
    """Reshape one published DE table into an AREE RNA-seq results file.

    `column_map` maps AREE column names to source column names. Any AREE column
    absent from the map is emitted empty rather than invented, so that
    downstream provenance shows exactly which statistics the source reported.

    Returns a report describing rows kept, rows dropped, and columns absent.
    """
    missing_source = [src for src in column_map.values() if src not in frame.columns]
    if missing_source:
        raise ValueError(f"source table lacks column(s): {missing_source}")

    out = pd.DataFrame()
    for aree_col in RNASEQ_COLUMNS:
        src = column_map.get(aree_col)
        out[aree_col] = frame[src] if src else pd.NA

    n_input = len(out)

    # Trailing annotation/footer rows carry no identifier or no effect size.
    out = out[out["gene_id"].notna() & out["log2FoldChange"].notna()]
    n_after_blank = len(out)

    # Published supplementary sheets are not always clean tables: a repeated
    # header row part-way down a sheet is common. Coerce the numeric columns and
    # drop anything that is not actually a number, counting it, rather than
    # letting a non-numeric string reach the harmonizer as an effect size.
    for col in ("log2FoldChange", "lfcSE", "pvalue", "padj"):
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out = out[out["log2FoldChange"].notna()]
    n_non_numeric = n_after_blank - len(out)

    out = out.reset_index(drop=True)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)

    return {
        "output_file": str(out_path.relative_to(REPO_ROOT))
        if out_path.is_relative_to(REPO_ROOT)
        else str(out_path),
        "rows_in_source": n_input,
        "rows_written": len(out),
        "rows_dropped_missing_id_or_effect": n_input - n_after_blank,
        "rows_dropped_non_numeric": n_non_numeric,
        "columns_reported_by_source": sorted(column_map),
        "columns_absent_from_source": sorted(set(RNASEQ_COLUMNS) - set(column_map)),
        "output_sha256": sha256_file(out_path),
    }


def write_intake_provenance(
    out_path: Path,
    *,
    study_id: str,
    source_description: str,
    source_url: str,
    source_license: str,
    source_file: Path,
    citation: str,
    conversions: list,
    notes: list,
    date_generated: str | None = None,
) -> None:
    """Record how a published artifact became AREE result files.

    `date_generated` records when these outputs were *derived*. Callers pass the
    existing date when re-running produces byte-identical outputs, so that a
    no-op re-run does not dirty the tree; it defaults to today otherwise.
    """
    doc = {
        "study_id": study_id,
        "date_generated": date_generated or date.today().isoformat(),
        "converter": "src/intake/supplementary_table.py",
        "source": {
            "description": source_description,
            "url": source_url,
            "license": source_license,
            "citation": citation,
            "local_copy": str(source_file.relative_to(REPO_ROOT))
            if source_file.is_relative_to(REPO_ROOT)
            else str(source_file),
            "sha256": sha256_file(source_file),
        },
        "conversions": conversions,
        "transformation_notes": notes,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(doc, indent=2) + "\n")
