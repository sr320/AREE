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
import warnings
from contextlib import contextmanager
from dataclasses import dataclass
from functools import lru_cache
from numbers import Integral, Real
from pathlib import Path

import pandas as pd

from common import (
    MAPPINGS_DIR,
    assembly_crosswalk_path,
    load_json,
    load_yaml,
    resolve_assembly,
)

# The demo crosswalk is the default so that `make demo` and the test suite run
# against self-contained synthetic data. A real study resolves against the
# crosswalk registered for its genome assembly in
# data/reference/genome_assemblies.yaml (`identifier_crosswalk`), so studies of
# different species in one run each use their own species' crosswalk. The
# AREE_CROSSWALK environment variable overrides that for local experiments, and
# only for a crosswalk built for the same taxon.
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
_CACHED_MTIMES: dict[str, int | None] = {}

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


def crosswalk_taxid(crosswalk_path) -> int | None:
    """The NCBI taxid a crosswalk was built for, from its provenance sidecar.

    None for the demo crosswalk and for any crosswalk without provenance.
    """
    crosswalk_path = Path(crosswalk_path)
    provenance = crosswalk_path.with_name(f"{crosswalk_path.stem}.provenance.json")
    if not provenance.exists():
        return None
    taxid = load_json(provenance).get("ncbi_taxid")
    return int(taxid) if taxid is not None else None


def crosswalk_route(crosswalk_path) -> str | None:
    """How a crosswalk was built, from its provenance: 'ncbi_gene' or 'gff_homology'.

    Crosswalks from build_crosswalk.py predate the field and are 'ncbi_gene';
    None for a crosswalk without provenance (the demo).
    """
    crosswalk_path = Path(crosswalk_path)
    provenance = crosswalk_path.with_name(f"{crosswalk_path.stem}.provenance.json")
    if not provenance.exists():
        return None
    return load_json(provenance).get("route", "ncbi_gene")


def crosswalk_for_study(study: dict) -> Path:
    """The crosswalk a study must resolve against, decided by the study itself.

    * A simulated study resolves against the synthetic demo crosswalk. Harmonizing
      the demo studies against a real crosswalk would silently attach real NCBI
      GeneIDs to fabricated effect sizes, because the demo LOC numbers collide
      with real genes.
    * A real study resolves against the crosswalk registered for its genome
      assembly (`identifier_crosswalk` in data/reference/genome_assemblies.yaml).
      The assembly, not the study's species, decides: identifiers come from the
      annotation the reads were mapped to.
    * $AREE_CROSSWALK overrides the registered crosswalk only if its provenance
      shows it was built for the assembly's taxon. A process-wide setting that
      could hand an oyster crosswalk to a mussel study is exactly what this
      function exists to prevent.

    The crosswalk last selected with `set_crosswalk_path()` is deliberately not
    consulted. It used to be, which meant the first real study harmonized in a
    process decided the crosswalk for every later one.

    A real study with no usable crosswalk is an error, never a silent fallback.
    """
    if study.get("simulated", False):
        return DEMO_CROSSWALK_PATH

    study_id = study.get("study_id", "<unknown>")
    assembly = resolve_assembly(study.get("genome_assembly"))
    if assembly is None:
        raise ValueError(
            f"study {study_id!r}: genome_assembly {study.get('genome_assembly')!r} is not "
            "in data/reference/genome_assemblies.yaml, so AREE cannot tell which "
            "crosswalk its identifiers belong to. Run `aree validate-study` first."
        )
    taxid = assembly.get("ncbi_taxid")
    registered = assembly_crosswalk_path(assembly)

    env = os.environ.get(CROSSWALK_ENV_VAR)
    if env:
        override = Path(env)
        override_taxid = crosswalk_taxid(override)
        if override_taxid is None:
            raise ValueError(
                f"${CROSSWALK_ENV_VAR}={env} has no provenance sidecar recording its taxid, "
                f"so it cannot be checked against study {study_id!r} (taxid {taxid}). "
                f"Unset ${CROSSWALK_ENV_VAR} to use the crosswalk registered for "
                f"{assembly['assembly_id']!r}."
            )
        if override_taxid != taxid:
            raise ValueError(
                f"${CROSSWALK_ENV_VAR}={env} was built for taxid {override_taxid}, but "
                f"study {study_id!r} is on {assembly['assembly_id']!r} (taxid {taxid}). "
                f"Unset ${CROSSWALK_ENV_VAR}; each study's crosswalk is read from "
                "data/reference/genome_assemblies.yaml."
            )
        if registered is not None and override.resolve() != registered.resolve():
            warnings.warn(
                f"${CROSSWALK_ENV_VAR} overrides the crosswalk registered for "
                f"{assembly['assembly_id']!r} ({registered}) with {override}.",
                stacklevel=2,
            )
        return override

    if registered is None:
        raise FileNotFoundError(
            f"study {study_id!r} is not simulated, and no identifier crosswalk is "
            f"registered for its assembly {assembly['assembly_id']!r} "
            f"({assembly.get('species')}, taxid {taxid}). Build one with\n"
            f"  aree build-crosswalk --taxid {taxid} --slug <short_name> "
            f"--organism \"{assembly.get('species')}\"\n"
            "and set `identifier_crosswalk` for this assembly in "
            "data/reference/genome_assemblies.yaml."
        )
    if not registered.exists():
        raise FileNotFoundError(
            f"study {study_id!r}: the crosswalk registered for {assembly['assembly_id']!r} "
            f"is missing: {registered}. Rebuild it with "
            f"`aree build-crosswalk --taxid {taxid}`."
        )
    return registered


def _clear_crosswalk_caches() -> None:
    _load_crosswalk.cache_clear()
    _column_index.cache_clear()
    _uniprot_index.cache_clear()
    _legacy_locus_index.cache_clear()
    _multi_id_index.cache_clear()
    _retired_map.cache_clear()


def set_crosswalk_path(path) -> None:
    """Select the crosswalk to resolve against, clearing cached lookups if the file changed.

    Pass None to fall back to $AREE_CROSSWALK / the demo crosswalk. This sets the
    crosswalk that `resolve_identifier` uses; it does not influence which
    crosswalk `crosswalk_for_study` chooses for the next study.
    """
    global _CROSSWALK_OVERRIDE
    new_path = Path(path) if path is not None else None
    _CROSSWALK_OVERRIDE = new_path
    if new_path is None:
        _CACHED_MTIMES.clear()
        _clear_crosswalk_caches()
        return
    key = str(new_path)
    mtime = new_path.stat().st_mtime_ns if new_path.exists() else None
    if key in _CACHED_MTIMES and _CACHED_MTIMES[key] != mtime:
        _clear_crosswalk_caches()
    _CACHED_MTIMES[key] = mtime


@contextmanager
def use_crosswalk(path):
    """Temporarily select a crosswalk for the duration of a block, restoring prior state."""
    global _CROSSWALK_OVERRIDE
    prev = _CROSSWALK_OVERRIDE
    set_crosswalk_path(path)
    try:
        yield
    finally:
        _CROSSWALK_OVERRIDE = prev


@lru_cache(maxsize=4)
def _load_crosswalk(path_str: str) -> pd.DataFrame:
    path = Path(path_str)
    if not path.exists():
        raise FileNotFoundError(
            f"Identifier crosswalk not found: {path}. Build one with "
            "`aree build-crosswalk`, or unset $AREE_CROSSWALK."
        )
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str).fillna("")
    if "gene_symbol" not in df.columns or not ({"ncbi_gene_id", "gene_id"} & set(df.columns)):
        raise ValueError(f"crosswalk {path} needs gene_symbol and an anchor column "
                         "(ncbi_gene_id, or gene_id for a GFF-route crosswalk)")
    if "orthogroup_id" not in df.columns:
        df["orthogroup_id"] = ""
    # The standardized identifier: NCBI GeneID for the NCBI route; the authors'
    # gene ID for the GFF route (mappings/gff_crosswalk.py), whose genes have no
    # GeneID. Everything downstream reads `_anchor`, never a route's own column.
    df["_anchor"] = df[anchor_column(df)]
    return df


def anchor_column(df: pd.DataFrame) -> str:
    """`gene_id` for a GFF-route crosswalk, else `ncbi_gene_id`."""
    return "gene_id" if "gene_id" in df.columns else "ncbi_gene_id"


# Multi-valued identifier columns of a GFF-route crosswalk: the authors'
# transcript and protein IDs, each naming exactly one gene in that annotation.
MULTI_ID_COLUMNS = ("transcript_ids", "protein_ids")


@lru_cache(maxsize=8)
def _multi_id_index(path_str: str, column: str) -> dict:
    """Map each ';'-delimited ID in `column` to its row index (-1 if two genes carry it)."""
    df = _load_crosswalk(path_str)
    index: dict = {}
    if column not in df.columns:
        return index
    for pos, ids in enumerate(df[column]):
        for value in (v for v in ids.split(";") if v):
            index[value] = -1 if value in index and index[value] != pos else pos
    return index


@lru_cache(maxsize=16)
def _column_index(path_str: str, column: str, lower: bool = False) -> dict[str, tuple[str, str | None, int]]:
    """Map each value in `column` to (first_anchor, first_orthogroup, match_count) for O(1) lookup.

    Preserves exact first-row semantics of `cw[cw[column] == candidate]`, including
    ';'-delimited multi-valued entries in `ensembl_gene_id` and skipping empty
    `ncbi_gene_id` values.
    """
    df = _load_crosswalk(path_str)
    index: dict[str, tuple[str, str | None, int]] = {}
    if column not in df.columns:
        return index
    split_semi = column == "ensembl_gene_id"
    skip_empty = column == "ncbi_gene_id"
    anchors = df["_anchor"].tolist()
    ogs = df["orthogroup_id"].tolist()
    vals = df[column].tolist()
    for val, anchor, og in zip(vals, anchors, ogs):
        if skip_empty and not val:
            continue
        og_val = og or None
        if split_semi:
            # Match both the full cell and each ';'-delimited token (deduped per row).
            tokens = {val} if val else set()
            tokens.update(t for t in val.split(";") if t)
            for tok in tokens:
                key = tok.lower() if lower else tok
                prev = index.get(key)
                index[key] = (anchor, og_val, 1) if prev is None else (prev[0], prev[1], prev[2] + 1)
        else:
            key = val.lower() if lower else val
            prev = index.get(key)
            index[key] = (anchor, og_val, 1) if prev is None else (prev[0], prev[1], prev[2] + 1)
    return index


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


@lru_cache(maxsize=4)
def _legacy_locus_index(path_str: str) -> dict:
    """Map each legacy locus tag (NCBI LocusTag, e.g. CGI_10001851) to its row index.

    Crosswalks built before builder v1.1.0 have no `legacy_locus_tags` column;
    they simply resolve nothing through this route. A tag carried by more than
    one gene maps to -1 and is treated as ambiguous.
    """
    df = _load_crosswalk(path_str)
    index: dict = {}
    if "legacy_locus_tags" not in df.columns:
        return index
    for pos, tags in enumerate(df["legacy_locus_tags"]):
        for tag in (t for t in tags.split(";") if t):
            index[tag] = -1 if tag in index and index[tag] != pos else pos
    return index


def crosswalk_annotation_release(crosswalk_path: Path) -> str | None:
    """The reference annotation the active crosswalk was built from, if recorded.

    This is the annotation every `feature_id_standardized` actually refers to,
    which is frequently *not* the assembly the source study used. Surfacing it
    per evidence record is what makes that crossing visible to a reader of the
    evidence table rather than something buried in a provenance sidecar.

    Returns None for the demo crosswalk, which is synthetic and has no
    provenance record.
    """
    crosswalk_path = Path(crosswalk_path)
    provenance = crosswalk_path.with_name(f"{crosswalk_path.stem}.provenance.json")
    if not provenance.exists():
        return None
    try:
        return load_json(provenance).get("ncbi_reference_annotation_at_build")
    except (ValueError, OSError):
        return None


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


# GFF/GTF-derived feature prefixes that published result tables commonly carry.
_FEATURE_PREFIXES = ("gene-", "rna-", "cds-", "exon-", "id-")


def stringify_identifier(raw_id) -> str:
    """Return a stable string representation for feature identifiers.

    Pandas row iteration may coerce an integer GeneID column to floats when the
    same row also contains floating-point statistics. A real NCBI GeneID like
    ``105337716`` can therefore arrive here as ``105337716.0`` even though the
    source table stored an integer identifier. Preserve the integral value in
    that case so exact crosswalk lookups keep working.
    """
    if isinstance(raw_id, Integral):
        return str(raw_id)
    if isinstance(raw_id, Real) and float(raw_id).is_integer():
        return str(int(raw_id))
    return str(raw_id).strip()


def identifier_candidates(raw_id: str) -> list:
    """Candidate tokens to look up for a raw identifier, most literal first.

    Published result tables rarely carry a bare identifier. Tables built from an
    NCBI GFF routinely report things like ``gene-LOC105320749|LOC105320749`` —
    a feature-type prefix plus a pipe-joined symbol. The decoration is
    deterministic and the stable identifier is literally embedded, so stripping
    it is lossless rather than a guess.

    The verbatim value is always tried first, and the caller keeps it as
    `feature_id_original`; this only affects what is looked up.
    """
    raw = stringify_identifier(raw_id)
    candidates = [raw]
    for token in raw.split("|"):
        token = token.strip()
        if not token:
            continue
        for prefix in _FEATURE_PREFIXES:
            if token.lower().startswith(prefix):
                token = token[len(prefix):]
                break
        if token and token not in candidates:
            candidates.append(token)
    return candidates


@lru_cache(maxsize=1)
def _ambiguous_map() -> dict:
    doc = load_yaml(AMBIGUOUS_MAP_PATH)
    return {entry["source_symbol"]: entry for entry in doc["entries"]}


def resolve_identifier(raw_id: str, id_type: str) -> ResolvedIdentifier:
    """Resolve an identifier against the active crosswalk, then attach its orthogroup.

    Real crosswalks ship `orthogroup_id` empty (it is not NCBI data). When an
    orthogroup run is active (data/reference/orthogroups/active.yaml), the
    resolved GeneID's orthogroup is taken from that run instead. The synthetic
    demo crosswalk keeps its own illustrative orthogroups and is never looked
    up in a real run.
    """
    resolved = _resolve_identifier(raw_id, id_type)
    # The active orthogroup run is keyed by NCBI GeneID, so only a GeneID-anchored
    # crosswalk may consult it: an author gene ID could collide with a GeneID.
    if (resolved.orthogroup_id is None and resolved.feature_id_standardized is not None
            and active_crosswalk_path() != DEMO_CROSSWALK_PATH
            and anchor_column(_crosswalk()) == "ncbi_gene_id"):
        from mappings.orthogroups import gene_orthogroup

        og = gene_orthogroup(resolved.feature_id_standardized)
        if og is not None:
            return ResolvedIdentifier(resolved.feature_id_standardized, resolved.mapping_confidence, og)
    return resolved


def _resolve_identifier(raw_id: str, id_type: str) -> ResolvedIdentifier:
    """Resolve a raw feature identifier to a standardized id + mapping confidence.

    id_type must be one of: ncbi_gene_id, ensembl_gene_id, uniprot_accession,
    locus_id, legacy_locus_tag, gene_symbol.
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
                row["_anchor"], "many_to_one_ortholog", row["orthogroup_id"] or None
            )
        row = cw.iloc[pos]
        # A GFF-route crosswalk reached this accession through our own homology
        # search, not a curated cross-reference: 'inferred', not 'exact'.
        confidence = "inferred" if row.get("uniprot_evidence", "") == "diamond_blastp" else "exact"
        return ResolvedIdentifier(row["_anchor"], confidence, row["orthogroup_id"] or None)

    if id_type == "legacy_locus_tag":
        # A legacy gene model (e.g. oyster_v9 CGI_ IDs) found in NCBI's LocusTag
        # for a current gene. 'inferred', not 'exact': the link crosses an
        # annotation version. A protein group naming several models resolves
        # only if every model points at the same current gene.
        index = _legacy_locus_index(str(active_crosswalk_path()))
        tags = [t.strip() for t in str(raw_id).replace(";", ",").split(",") if t.strip()]
        positions = {index.get(t) for t in tags}
        if not tags or None in positions or -1 in positions or len(positions) != 1:
            return ResolvedIdentifier(None, "unresolved", None)
        row = cw.iloc[positions.pop()]
        return ResolvedIdentifier(row["_anchor"], "inferred", row["orthogroup_id"] or None)

    if id_type in STABLE_ID_COLUMNS:
        # `ncbi_gene_id` and `locus_id` are two spellings of the same thing:
        # NCBI's LOC<GeneID> convention means "LOC105320749" and "105320749"
        # denote one gene. Published tables use whichever their annotation
        # pipeline emitted, so both columns are searched for either id_type
        # rather than trusting the caller's guess from a column name.
        if id_type in ("ncbi_gene_id", "locus_id"):
            columns = ["ncbi_gene_id", "locus_id"]
            if "gene_id" in cw.columns:
                columns = ["gene_id", "locus_id", "ncbi_gene_id"]
        else:
            columns = [STABLE_ID_COLUMNS[id_type]]
        columns = [c for c in columns if c in cw.columns]
        if not columns:
            return ResolvedIdentifier(None, "unresolved", None)

        cw_path_str = str(active_crosswalk_path())
        candidates = identifier_candidates(raw_id)
        for candidate in candidates:
            for column in columns:
                hit = _column_index(cw_path_str, column).get(candidate)
                if hit is not None:
                    anchor, og, count = hit
                    confidence = "exact" if count == 1 else "many_to_one_ortholog"
                    return ResolvedIdentifier(anchor, confidence, og)
            # GFF route: the authors' transcript or protein ID names one gene.
            if id_type in ("ncbi_gene_id", "locus_id"):
                for column in MULTI_ID_COLUMNS:
                    pos = _multi_id_index(cw_path_str, column).get(candidate)
                    if pos is not None and pos >= 0:
                        row = cw.iloc[pos]
                        return ResolvedIdentifier(row["_anchor"], "exact", row["orthogroup_id"] or None)

        # Nothing in the current annotation. The identifier may be a GeneID that
        # NCBI has since retired — common for studies run on an older annotation.
        if id_type in ("ncbi_gene_id", "locus_id"):
            retired = _retired_map(cw_path_str)
            ncbi_idx = _column_index(cw_path_str, "ncbi_gene_id")
            for candidate in candidates:
                probe = candidate[3:] if candidate.startswith("LOC") else candidate
                current = retired.get(probe)
                if not current:
                    continue
                hit = ncbi_idx.get(current)
                if hit is not None and hit[2] == 1:
                    # 'inferred', not 'exact': the remap crosses an annotation version.
                    return ResolvedIdentifier(hit[0], "inferred", hit[1])
        return ResolvedIdentifier(None, "unresolved", None)

    if id_type == "gene_symbol":
        cw_path_str = str(active_crosswalk_path())
        hit = _column_index(cw_path_str, "gene_symbol", lower=True).get(str(raw_id).lower())
        if hit is not None:
            anchor, og, count = hit
            # Symbol-based match to a stable reference id, but symbols are not
            # guaranteed unique/stable, so this is "inferred" rather than "exact".
            confidence = "inferred" if count == 1 else "many_to_one_ortholog"
            return ResolvedIdentifier(anchor, confidence, og)

        ambiguous = _ambiguous_map().get(raw_id)
        if ambiguous:
            standardized = ambiguous["standardized_id"]
            og_hit = _column_index(cw_path_str, "_anchor").get(standardized)
            orthogroup = og_hit[1] if og_hit is not None else None
            return ResolvedIdentifier(standardized, ambiguous["mapping_confidence"], orthogroup)

        return ResolvedIdentifier(None, "unresolved", None)

    raise ValueError(f"Unknown id_type: {id_type!r}")
