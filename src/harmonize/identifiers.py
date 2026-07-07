"""Identifier harmonization.

Precedence hierarchy (highest to lowest trust), per docs/identifier_mapping.md:
    NCBI Gene ID -> Ensembl gene ID -> UniProt accession -> reference locus ID
    -> gene symbol -> orthogroup

The original identifier reported by the source study is always preserved by
the caller (feature_id_original); this module only computes the standardized
identifier and a mapping_confidence.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import pandas as pd

from common import MAPPINGS_DIR, load_yaml

CROSSWALK_PATH = MAPPINGS_DIR / "gene_id_crosswalk.tsv"
AMBIGUOUS_MAP_PATH = MAPPINGS_DIR / "ambiguous_symbol_map.yaml"

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


@lru_cache(maxsize=1)
def _crosswalk() -> pd.DataFrame:
    return pd.read_csv(CROSSWALK_PATH, sep="\t", comment="#", dtype=str)


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

    if id_type in STABLE_ID_COLUMNS:
        column = STABLE_ID_COLUMNS[id_type]
        match = cw[cw[column] == raw_id]
        if len(match) == 1:
            row = match.iloc[0]
            return ResolvedIdentifier(row["ncbi_gene_id"], "exact", row["orthogroup_id"])
        if len(match) > 1:
            row = match.iloc[0]
            return ResolvedIdentifier(row["ncbi_gene_id"], "many_to_one_ortholog", row["orthogroup_id"])
        return ResolvedIdentifier(None, "unresolved", None)

    if id_type == "gene_symbol":
        match = cw[cw["gene_symbol"].str.lower() == str(raw_id).lower()]
        if len(match) == 1:
            row = match.iloc[0]
            # Symbol-based match to a stable reference id, but symbols are not
            # guaranteed unique/stable, so this is "inferred" rather than "exact".
            return ResolvedIdentifier(row["ncbi_gene_id"], "inferred", row["orthogroup_id"])
        if len(match) > 1:
            row = match.iloc[0]
            return ResolvedIdentifier(row["ncbi_gene_id"], "many_to_one_ortholog", row["orthogroup_id"])

        ambiguous = _ambiguous_map().get(raw_id)
        if ambiguous:
            standardized = ambiguous["standardized_id"]
            og_row = cw[cw["ncbi_gene_id"] == standardized]
            orthogroup = og_row.iloc[0]["orthogroup_id"] if len(og_row) else None
            return ResolvedIdentifier(standardized, ambiguous["mapping_confidence"], orthogroup)

        return ResolvedIdentifier(None, "unresolved", None)

    raise ValueError(f"Unknown id_type: {id_type!r}")
