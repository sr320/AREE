"""Export AREE's queryable tables as Parquet for the site's in-browser explorer.

`aree export-site-data` writes docs/explore/data/:

* evidence.parquet   — every harmonized evidence record (the EVIDENCE_COLUMNS
                       subset), sorted by standardized feature id so a gene
                       lookup touches only the row groups holding that gene;
* features.parquet   — one row per distinct identifier combination (standardized,
                       original, orthogroup), the small index a lookup searches first;
* candidates.parquet — every ranked candidate, with gene symbol and description
                       for real genes (from the built crosswalks);
* studies.parquet    — one row per registered study;
* meta.json          — build date, row counts, real/simulated split, and the
                       interpretation note shown on the page.

docs/explore.qmd queries these with DuckDB-WASM in the visitor's browser. DuckDB
reads Parquet over HTTP range requests, fetching only the columns and row groups
a query needs, so the full evidence table never has to be downloaded.

Requires pyarrow (the `site` extra): pip install -e ".[site]"
"""
from __future__ import annotations

import datetime as _dt
import json
from pathlib import Path

import pandas as pd

import query
from common import REPO_ROOT, load_vocab

DEFAULT_OUT_DIR = REPO_ROOT / "docs" / "explore" / "data"

# Small row groups keep a single-gene lookup to a few hundred KB of transfer.
ROW_GROUP_SIZE = 20_000

# Columns written as real booleans, so `WHERE is_predictive_baseline` works in SQL.
BOOLEAN_COLUMNS = {"simulated", "is_predictive_baseline"}


def _normalize(df: pd.DataFrame) -> pd.DataFrame:
    """One consistent Parquet type per column: bool, float, or nullable text."""
    df = df.copy()
    for col in df.columns:
        if col in BOOLEAN_COLUMNS:
            df[col] = df[col].astype(str).str.strip().str.lower().isin({"true", "1", "yes"})
        elif col in query.NUMERIC_COLUMNS:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        elif df[col].dtype == object:
            df[col] = df[col].map(lambda v: None if v is None or (isinstance(v, float) and pd.isna(v)) else str(v))
    return df


def _with_gene_names(candidates: pd.DataFrame) -> pd.DataFrame:
    """Add gene_symbol and gene_description for gene candidates whose ID is an NCBI
    GeneID in a built crosswalk; blank otherwise (demo IDs, proteins, metabolites)."""
    from reporting.dashboard import gene_annotation_frame

    names = gene_annotation_frame()
    ids = candidates["feature_id_standardized"].astype(str)
    is_gene = candidates["feature_type"].astype(str) == "gene"
    by_id = names.set_index("ncbi_gene_id")
    out = candidates.copy()
    for col in ("gene_symbol", "gene_description"):
        out[col] = ids.map(by_id[col]).where(is_gene) if len(by_id) else None
    return out


def _write(df: pd.DataFrame, path: Path) -> None:
    try:
        df.to_parquet(path, index=False, compression="zstd", row_group_size=ROW_GROUP_SIZE)
    except ImportError as exc:
        raise SystemExit('Exporting site data needs pyarrow: pip install -e ".[site]"') from exc


def _split(df: pd.DataFrame) -> dict:
    if df.empty or "simulated" not in df.columns:
        return {"total": len(df), "real": 0, "simulated": 0}
    sim = df["simulated"].astype(bool)
    return {"total": len(df), "real": int((~sim).sum()), "simulated": int(sim.sum())}


def export_site_data(out_dir: Path = DEFAULT_OUT_DIR,
                     sources: query.Sources = query.DEFAULT_SOURCES,
                     include_simulated: bool = False) -> dict:
    """Write the explorer's Parquet tables and meta.json; return the meta document.

    Simulated demo studies and their evidence are left out unless
    include_simulated is set: they are pipeline fixtures, not findings.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    evidence = sources.evidence()
    if not evidence.empty:
        evidence = evidence[[c for c in query.EVIDENCE_COLUMNS if c in evidence.columns]]
        # Sorted by standardized id, then original id (which orders the unresolved
        # rows at the end), so exact-id lookups prune to a few row groups.
        evidence = _normalize(evidence).sort_values(
            ["feature_id_standardized", "feature_id_original", "study_id"], na_position="last")
    else:
        evidence = pd.DataFrame(columns=query.EVIDENCE_COLUMNS)

    candidates = sources.candidates()
    if not candidates.empty:
        candidates = candidates[[c for c in query.CANDIDATE_COLUMNS if c in candidates.columns]].copy()
        candidates["card_file"] = candidates["card_file"].map(
            lambda p: Path(p).name if isinstance(p, str) and p else None)
        candidates = _normalize(_with_gene_names(candidates))
    else:
        candidates = pd.DataFrame(columns=query.CANDIDATE_COLUMNS)

    studies = _normalize(query.study_table(sources))

    evidence = query._filter_simulated(evidence, include_simulated)
    candidates = query._filter_simulated(candidates, include_simulated)
    studies = query._filter_simulated(studies, include_simulated)
    tables = {"evidence": evidence, "features": _feature_index(evidence),
              "candidates": candidates, "studies": studies}
    for name, df in tables.items():
        _write(df.reset_index(drop=True), out_dir / f"{name}.parquet")

    meta = {
        "generated": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "aree_version": _aree_version(),
        "tables": {name: {"file": f"{name}.parquet", "columns": list(df.columns), **_split(df)}
                   for name, df in tables.items()},
        "facets": {
            "evidence": _facets(evidence, ["phenotype", "stressor", "species", "feature_type", "tissue",
                                           "life_stage", "mapping_confidence"]),
            "candidates": _facets(candidates, ["phenotype", "tier", "evidence_class", "contrast_type",
                                               "feature_type"]),
        },
        "labels": {
            "phenotype": {tid: t.get("label", tid) for tid, t in load_vocab("phenotype_ontology").items()},
            "stressor": {tid: t.get("label", tid) for tid, t in load_vocab("stressor_ontology").items()},
        },
        "artifacts": query.artifacts(sources),
        "note": query.INTERPRETATION_NOTE,
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    return meta


FEATURE_INDEX_COLUMNS = ["feature_id_standardized", "feature_id_original", "orthogroup_id",
                         "feature_type", "species", "simulated"]


def _feature_index(evidence: pd.DataFrame) -> pd.DataFrame:
    """One row per distinct identifier combination, with its record count.

    The page searches this small table (substring matches included) to resolve
    what the user typed, then fetches evidence by exact id, so a lookup never
    scans the full evidence file.
    """
    cols = [c for c in FEATURE_INDEX_COLUMNS if c in evidence.columns]
    if evidence.empty:
        return pd.DataFrame(columns=[*FEATURE_INDEX_COLUMNS, "n_records"])
    return (evidence[cols].astype({c: object for c in cols if c != "simulated"})
            .fillna({c: "" for c in cols if c != "simulated"})
            .groupby(cols, dropna=False).size().rename("n_records").reset_index()
            .replace({"": None})
            .sort_values(["feature_id_standardized", "feature_id_original"], na_position="last"))


def _facets(df: pd.DataFrame, columns: list[str]) -> dict:
    """Distinct values per filter column, for the page's dropdowns."""
    return {c: sorted(df[c].dropna().astype(str).unique().tolist()) for c in columns if c in df.columns}


def _aree_version() -> str:
    from aree import __version__
    return __version__
