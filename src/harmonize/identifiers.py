"""Identifier harmonization.

Precedence hierarchy (highest to lowest trust), per docs/identifier_mapping.md:
    NCBI Gene ID -> Ensembl gene ID -> UniProt accession -> reference locus ID
    -> gene symbol -> orthogroup

The original identifier reported by the source study is always preserved by
the caller (feature_id_original); this module only computes the standardized
identifier and a mapping_confidence.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import pandas as pd

from common import MAPPINGS_DIR, load_yaml

# The demo crosswalk is the default so that `make demo` and the test suite run
# against self-contained synthetic data. Curating a real study requires pointing
# at a real crosswalk built by src/mappings/build_crosswalk.py, either via the
# AREE_CROSSWALK environment variable or `set_crosswalk_path()`.
#
# These two files are deliberately NOT merged: the demo LOC numbers collide with
# real NCBI GeneIDs while carrying different (fabricated) gene identities, so a
# union would silently attach demo biology to real accessions.
DEMO_CROSSWALK_PATH = MAPPINGS_DIR / "gene_id_crosswalk.tsv"
CROSSWALK_ENV_VAR = "AREE_CROSSWALK"
AMBIGUOUS_MAP_PATH = MAPPINGS_DIR / "ambiguous_symbol_map.yaml"

# Backwards-compatible alias; prefer active_crosswalk_path().
CROSSWALK_PATH = DEMO_CROSSWALK_PATH

_CROSSWALK_OVERRIDE: Path | None = None

STABLE_ID_COLUMNS = {
    "ncbi_gene_id": "ncbi_gene_id",
    "ensembl_gene_id": "ensembl_gene_id",
    "uniprot_accession": "uniprot_accession",
    "locus_id": "locus_id",
}


@dataclass(frozen=True)
class ResolvedIdentifier:
    feature_id_standardized: str | None
    mapping_confidence: str
    orthogroup_id: str | None


def active_crosswalk_path() -> Path:
    """Path of the crosswalk currently in force.

    Precedence: explicit `set_crosswalk_path()` > $AREE_CROSSWALK > demo crosswalk.
    """
    if _CROSSWALK_OVERRIDE is not None:
        return _CROSSWALK_OVERRIDE
    env = os.environ.get(CROSSWALK_ENV_VAR)
    if env:
        return Path(env)
    return DEMO_CROSSWALK_PATH


def set_crosswalk_path(path) -> None:
    """Select the crosswalk to resolve against, clearing cached lookups.

    Pass None to fall back to $AREE_CROSSWALK / the demo crosswalk.
    """
    global _CROSSWALK_OVERRIDE
    _CROSSWALK_OVERRIDE = Path(path) if path is not None else None
    _load_crosswalk.cache_clear()
    _uniprot_index.cache_clear()
    _retired_map.cache_clear()


@lru_cache(maxsize=4)
def _load_crosswalk(path_str: str) -> pd.DataFrame:
    path = Path(path_str)
    if not path.exists():
        raise FileNotFoundError(
            f"Identifier crosswalk not found: {path}. Build one with "
            "`aree build-crosswalk`, or unset $AREE_CROSSWALK to use the demo crosswalk."
        )
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str).fillna("")
    missing = {"ncbi_gene_id", "gene_symbol"} - set(df.columns)
    if missing:
        raise ValueError(f"crosswalk {path} is missing required column(s): {sorted(missing)}")
    if "orthogroup_id" not in df.columns:
        df["orthogroup_id"] = ""
    return df


def _crosswalk() -> pd.DataFrame:
    return _load_crosswalk(str(active_crosswalk_path()))


@lru_cache(maxsize=4)
def _uniprot_index(path_str: str) -> dict:
    """Map every UniProt accession for a gene (not just the representative) to its row index.

    Real genes carry multiple TrEMBL accessions, so the crosswalk keeps the full
    set in `uniprot_accessions_all` while `uniprot_accession` holds one
    representative. Without this index, a study reporting a non-representative
    accession would resolve as `unresolved` even though the link is exact.
    """
    df = _load_crosswalk(path_str)
    index: dict = {}
    collisions = set()
    for pos, row in enumerate(df.itertuples(index=False)):
        accessions = set()
        primary = getattr(row, "uniprot_accession", "") or ""
        if primary:
            accessions.add(primary)
        all_field = getattr(row, "uniprot_accessions_all", "") or ""
        accessions.update(a for a in all_field.split(";") if a)
        for acc in accessions:
            if acc in index and index[acc] != pos:
                collisions.add(acc)
            index.setdefault(acc, pos)
    # An accession attributed to more than one gene cannot be an exact match.
    for acc in collisions:
        index[acc] = -1
    return index


def retired_table_path(crosswalk_path: Path) -> Path | None:
    """Sidecar path for a crosswalk, by convention.

    `<slug>_gene_id_crosswalk.tsv` -> `<slug>_retired_gene_ids.tsv`.

    Returns None when the crosswalk does not follow the slug convention (the demo
    crosswalk is a bare `gene_id_crosswalk.tsv`), so that the naive substitution
    cannot resolve back onto the crosswalk file itself.
    """
    name = crosswalk_path.name
    sidecar = name.replace("_gene_id_crosswalk.tsv", "_retired_gene_ids.tsv")
    if sidecar == name:
        return None
    return crosswalk_path.with_name(sidecar)


@lru_cache(maxsize=4)
def _retired_map(path_str: str) -> dict:
    """Map discontinued NCBI GeneIDs to the current GeneID that replaced them.

    Studies analysed against an older annotation report GeneIDs that NCBI has
    since retired. Without this, such identifiers resolve as `unresolved` even
    though NCBI records an authoritative replacement. GeneIDs discontinued with
    no replacement are absent from the sidecar and stay unresolved by design.

    Returns an empty mapping when no sidecar exists (e.g. the demo crosswalk).
    """
    path = retired_table_path(Path(path_str))
    if path is None or not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str).fillna("")
    return dict(zip(df["discontinued_gene_id"], df["current_gene_id"]))


@lru_cache(maxsize=1)
def _ambiguous_map() -> dict:
    doc = load_yaml(AMBIGUOUS_MAP_PATH)
    return {entry["source_symbol"]: entry for entry in doc["entries"]}


def resolve_identifier(raw_id: str, id_type: str) -> ResolvedIdentifier:
    """Resolve a raw feature identifier to a standardized id + mapping confidence.

    id_type must be one of: ncbi_gene_id, ensembl_gene_id, uniprot_accession,
    locus_id, gene_symbol.
    """
    if raw_id is None or (isinstance(raw_id, float) and pd.isna(raw_id)):
        return ResolvedIdentifier(None, "unresolved", None)

    cw = _crosswalk()

    if id_type == "uniprot_accession":
        pos = _uniprot_index(str(active_crosswalk_path())).get(str(raw_id))
        if pos is None:
            return ResolvedIdentifier(None, "unresolved", None)
        if pos == -1:
            # Same accession attributed to multiple genes.
            match = cw[cw["uniprot_accessions_all"].str.contains(str(raw_id), regex=False)]
            row = match.iloc[0]
            return ResolvedIdentifier(
                row["ncbi_gene_id"], "many_to_one_ortholog", row["orthogroup_id"] or None
            )
        row = cw.iloc[pos]
        return ResolvedIdentifier(row["ncbi_gene_id"], "exact", row["orthogroup_id"] or None)

    if id_type in STABLE_ID_COLUMNS:
        column = STABLE_ID_COLUMNS[id_type]
        if column not in cw.columns:
            return ResolvedIdentifier(None, "unresolved", None)
        # Ensembl xrefs can be multi-valued (';'-delimited) in a real crosswalk.
        if column == "ensembl_gene_id":
            match = cw[
                (cw[column] == raw_id)
                | cw[column].str.split(";").apply(lambda ids: raw_id in ids)
            ]
        else:
            match = cw[cw[column] == raw_id]
        if len(match) == 1:
            row = match.iloc[0]
            return ResolvedIdentifier(row["ncbi_gene_id"], "exact", row["orthogroup_id"] or None)
        if len(match) > 1:
            row = match.iloc[0]
            return ResolvedIdentifier(row["ncbi_gene_id"], "many_to_one_ortholog", row["orthogroup_id"] or None)

        # Nothing in the current annotation. The identifier may be a GeneID that
        # NCBI has since retired — common for studies run on an older annotation.
        if id_type in ("ncbi_gene_id", "locus_id"):
            probe = str(raw_id)
            if id_type == "locus_id" and probe.startswith("LOC"):
                probe = probe[3:]
            current = _retired_map(str(active_crosswalk_path())).get(probe)
            if current:
                repl = cw[cw["ncbi_gene_id"] == current]
                if len(repl) == 1:
                    row = repl.iloc[0]
                    # 'inferred', not 'exact': the remap crosses an annotation version.
                    return ResolvedIdentifier(
                        row["ncbi_gene_id"], "inferred", row["orthogroup_id"] or None
                    )
        return ResolvedIdentifier(None, "unresolved", None)

    if id_type == "gene_symbol":
        match = cw[cw["gene_symbol"].str.lower() == str(raw_id).lower()]
        if len(match) == 1:
            row = match.iloc[0]
            # Symbol-based match to a stable reference id, but symbols are not
            # guaranteed unique/stable, so this is "inferred" rather than "exact".
            return ResolvedIdentifier(row["ncbi_gene_id"], "inferred", row["orthogroup_id"] or None)
        if len(match) > 1:
            row = match.iloc[0]
            return ResolvedIdentifier(row["ncbi_gene_id"], "many_to_one_ortholog", row["orthogroup_id"] or None)

        ambiguous = _ambiguous_map().get(raw_id)
        if ambiguous:
            standardized = ambiguous["standardized_id"]
            og_row = cw[cw["ncbi_gene_id"] == standardized]
            orthogroup = (og_row.iloc[0]["orthogroup_id"] or None) if len(og_row) else None
            return ResolvedIdentifier(standardized, ambiguous["mapping_confidence"], orthogroup)

        return ResolvedIdentifier(None, "unresolved", None)

    raise ValueError(f"Unknown id_type: {id_type!r}")
