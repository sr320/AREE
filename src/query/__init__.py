"""Read-only query layer over AREE's registry and pipeline outputs.

One set of functions answers "what does AREE hold?" questions — studies,
harmonized evidence, ranked candidates, evidence cards, vocabularies — for the
MCP server (`query.mcp_server`) and the static-site export (`query.site_export`).

It reads what the `aree` pipeline has already written; it never recomputes a
meta-analysis or writes anything. Every function returns plain JSON-serializable
dicts/lists, and every record keeps its provenance (study_id, evidence_id,
source_file, workflow_version, date_generated) and its `simulated` flag, so a
caller — human or language model — can always say where a number came from.

Evidence and candidates are queried with DuckDB, from either of two sources
(`Sources`): the pipeline's TSVs under reports/ (a curator's checkout), or the
Parquet tables `aree export-site-data` publishes (a server, or anyone who ran
`aree fetch-site-data`), which DuckDB reads in place without loading them.
Study records always come from the registry YAMLs.

Missing pipeline outputs are reported as missing (`artifacts`), never silently
treated as "no evidence".
"""
from __future__ import annotations

import math
import threading
from dataclasses import dataclass, field
from pathlib import Path

import duckdb
import pandas as pd

from common import (
    DEMO_STUDIES_DIR,
    REPO_ROOT,
    REPORTS_DIR,
    STUDIES_DIR,
    load_vocab,
    load_yaml,
    resolve_species,
)

# Stated with every answer that carries evidence: AREE reports associations and
# convergence, never validation or causation (CLAUDE.md, docs/interpreting_candidate_scores.md).
INTERPRETATION_NOTE = (
    "AREE evidence records are statistical associations from individual studies. "
    "A significant record or a high candidate score is not a validated biomarker and "
    "does not establish causation; emerging-tier candidates require replication. "
    "Records with simulated=true come from synthetic demo studies and are not real findings."
)

# Columns returned for an evidence record: identity, effect, context, and provenance.
EVIDENCE_COLUMNS = [
    "evidence_id", "study_id", "comparison_id", "contrast_type", "simulated", "species",
    "feature_id_original", "feature_id_standardized", "feature_type", "orthogroup_id",
    "mapping_confidence", "molecular_direction", "effect_size", "effect_size_type",
    "standard_error", "ci_lower", "ci_upper", "p_value", "adjusted_p_value", "sample_size",
    "tissue", "life_stage", "stressor", "phenotype", "phenotype_direction",
    "resilience_classification", "analysis_method", "quality_flags",
    "genome_assembly", "annotation_version", "source_file", "workflow_version", "date_generated",
]

CANDIDATE_COLUMNS = [
    "feature_id_standardized", "feature_type", "phenotype", "simulated", "species_taxid",
    "tier", "evidence_class", "contrast_type", "is_predictive_baseline", "context_replication",
    "score", "k_studies",
    "n_supporting_layers", "direction_consistency", "pooled_effect", "i_squared",
    "adjusted_p_value", "studies", "card_file",
]

STUDY_COLUMNS = [
    "study_id", "simulated", "species", "species_accepted", "assay_type", "analysis_mode",
    "tissues", "life_stages", "stressors", "phenotypes", "resilience_classifications",
    "contrast_types", "n_comparisons", "data_status", "raw_data_available", "processed_data_available",
    "qc_status", "analysis_status", "doi", "bioproject", "citation",
]

MAX_LIMIT = 500

# Read as text, then coerced: identifiers such as NCBI Gene IDs look numeric but must
# stay strings (the same rule harmonize.core applies to the evidence table).
NUMERIC_COLUMNS = {
    "effect_size", "standard_error", "ci_lower", "ci_upper", "p_value", "adjusted_p_value",
    "sample_size", "score", "k_studies", "n_supporting_layers", "direction_consistency",
    "pooled_effect", "i_squared", "species_taxid", "n_comparisons",
}

BOOLEAN_COLUMNS = {"simulated", "is_predictive_baseline"}

# Pairs of molecular directions that contradict each other: expression or abundance
# (up/down) and methylation (hyper/hypo).
OPPOSING_DIRECTIONS = (("up", "down"), ("hyper", "hypo"))

SIGNIFICANCE_THRESHOLD = 0.05


def _read_tsv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    for col in NUMERIC_COLUMNS & set(df.columns):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _sql_literal(path: Path) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def _attach(con: duckdb.DuckDBPyConnection, name: str, path: Path) -> None:
    """Expose one output as table `name`, with `_row` keeping the file's row order.

    Parquet is read in place (a view); a TSV is parsed exactly as before (text, then
    numeric columns coerced) and copied into DuckDB, so both sources answer alike.
    """
    if path.suffix == ".parquet":
        con.execute(f"CREATE VIEW {name} AS SELECT * EXCLUDE (file_row_number), "
                    f"file_row_number AS _row FROM read_parquet({_sql_literal(path)}, file_row_number = true)")
        return
    frame = _read_tsv(path)
    frame["_row"] = range(len(frame))
    con.register("_frame", frame)
    con.execute(f"CREATE TABLE {name} AS SELECT * FROM _frame")
    con.unregister("_frame")


def _stamp(path: Path):
    return path.stat().st_mtime_ns if path.exists() else None


@dataclass
class Sources:
    """Where the query layer reads from; tests point these at a temp directory.

    By default the pipeline's outputs under `reports_dir`. With `data_dir` set,
    the Parquet tables written by `aree export-site-data` instead (evidence,
    candidates and evidence cards), so answers need no pipeline run.
    """

    reports_dir: Path = REPORTS_DIR
    studies_dir: Path = STUDIES_DIR
    # Simulated demo studies live apart from real ones; `include_simulated`
    # decides whether answers show them, as for evidence and candidates.
    demo_studies_dir: Path = DEMO_STUDIES_DIR
    data_dir: Path | None = None
    _cache: dict = field(default_factory=dict, repr=False)
    _lock: threading.Lock = field(default_factory=threading.Lock, repr=False)
    _con: duckdb.DuckDBPyConnection | None = field(default=None, repr=False)
    _con_stamp: tuple | None = field(default=None, repr=False)

    @classmethod
    def from_site_data(cls, data_dir) -> Sources:
        return cls(data_dir=Path(data_dir))

    def study_paths(self) -> list[Path]:
        """Every study YAML, real then demo; templates (`_*`) skipped."""
        return [p for d in (self.studies_dir, self.demo_studies_dir)
                for p in sorted(d.glob("*.yaml")) if not p.name.startswith("_")]

    @property
    def evidence_table(self) -> Path:
        if self.data_dir is not None:
            return self.data_dir / "evidence.parquet"
        return self.reports_dir / "evidence" / "evidence_table.tsv"

    @property
    def candidates_table(self) -> Path:
        if self.data_dir is not None:
            return self.data_dir / "candidates.parquet"
        return self.reports_dir / "evidence_cards" / "candidates.tsv"

    @property
    def cards_dir(self) -> Path:
        return self.reports_dir / "evidence_cards"

    @property
    def cards_table(self) -> Path | None:
        return self.data_dir / "cards.parquet" if self.data_dir is not None else None

    def _read(self, path: Path) -> pd.DataFrame:
        if not path.exists():
            return pd.DataFrame()
        key = (str(path), path.stat().st_mtime_ns)
        if key not in self._cache:
            self._cache[key] = pd.read_parquet(path) if path.suffix == ".parquet" else _read_tsv(path)
        return self._cache[key]

    def evidence(self) -> pd.DataFrame:
        """The whole evidence table as a DataFrame (for export; queries use `db`)."""
        return self._read(self.evidence_table)

    def candidates(self) -> pd.DataFrame:
        return self._read(self.candidates_table)

    def db(self) -> duckdb.DuckDBPyConnection:
        """A cursor on an in-memory DuckDB holding `evidence` and `candidates` (each
        absent when its output is missing), rebuilt when an output file changes.
        One cursor per call, so concurrent callers (a threaded server) don't share one."""
        tables = (("evidence", self.evidence_table), ("candidates", self.candidates_table))
        stamp = tuple(_stamp(path) for _, path in tables)
        with self._lock:
            if self._con is None or self._con_stamp != stamp:
                con = duckdb.connect()
                for name, path in tables:
                    if path.exists():
                        _attach(con, name, path)
                self._con, self._con_stamp = con, stamp
            return self._con.cursor()


DEFAULT_SOURCES = Sources()


# --------------------------------------------------------------------------- helpers

def _clean(value):
    """NaN/inf/NA -> None, numpy scalars -> Python, so results serialize as JSON."""
    if value is None or value is pd.NA:
        return None
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _records(df: pd.DataFrame, columns: list[str]) -> list[dict]:
    cols = [c for c in columns if c in df.columns]
    return [{k: _clean(v) for k, v in row.items()} for row in df[cols].to_dict("records")]


def _limit(limit) -> int:
    return max(1, min(int(limit or 50), MAX_LIMIT))


def _as_bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def _with_booleans(df: pd.DataFrame) -> pd.DataFrame:
    """Flags as real booleans whichever source they came from (TSV text or Parquet bool)."""
    for col in BOOLEAN_COLUMNS & set(df.columns):
        df[col] = _as_bool_series(df[col])
    return df


def _contains(series: pd.Series, needle: str) -> pd.Series:
    return series.astype(str).str.lower().str.contains(needle.lower(), na=False, regex=False)


def _filter_simulated(df: pd.DataFrame, include_simulated: bool) -> pd.DataFrame:
    if include_simulated or "simulated" not in df.columns or df.empty:
        return df
    return df[~_as_bool_series(df["simulated"])]


def _conflicting(directions) -> bool:
    present = set(directions)
    return any(a in present and b in present for a, b in OPPOSING_DIRECTIONS)


def _repo_relative(path) -> str:
    try:
        return str(Path(path).resolve().relative_to(REPO_ROOT))
    except (ValueError, OSError):
        return str(path)


# SQL over columns that are text when read from a TSV and typed when read from Parquet.
def _text(column: str) -> str:
    return f"lower(CAST({column} AS VARCHAR))"


def _truthy(column: str) -> str:
    return f"coalesce(trim({_text(column)}) IN ('true', '1', 'yes'), false)"


def _number(column: str) -> str:
    return f"TRY_CAST({column} AS DOUBLE)"


class _Table:
    """One queryable table: its columns, and a WHERE clause built up with bound parameters."""

    def __init__(self, sources: Sources, name: str):
        self.name = name
        self.cur = sources.db()
        try:
            self.cur.execute(f"SELECT * FROM {name} LIMIT 0")
            self.columns = {d[0] for d in self.cur.description} - {"_row"}
            self.n_rows = self.cur.execute(f"SELECT count(*) FROM {name}").fetchone()[0]
        except duckdb.CatalogException:
            self.columns, self.n_rows = set(), 0
        self.clauses: list[str] = []
        self.params: list = []

    def where(self, clause: str, *params) -> _Table:
        self.clauses.append(clause)
        self.params.extend(params)
        return self

    def equals(self, column: str, value) -> _Table:
        """Case-insensitive equality; a filter on a column this table lacks is ignored."""
        if column not in self.columns:
            return self
        return self.where(f"{_text(column)} = lower(?)", str(value))

    def contains(self, column: str, value) -> _Table:
        return self.where(f"contains({_text(column)}, lower(?))", str(value))

    def feature(self, value: str, exact: bool) -> _Table:
        """Match a standardized or original identifier, or an orthogroup."""
        cols = [c for c in ("feature_id_standardized", "feature_id_original", "orthogroup_id")
                if c in self.columns]
        test = "{} = lower(?)" if exact else "contains({}, lower(?))"
        return self.where("(" + " OR ".join(test.format(_text(c)) for c in cols) + ")",
                          *[value] * len(cols))

    def real_only(self, include_simulated: bool) -> _Table:
        if include_simulated or "simulated" not in self.columns:
            return self
        return self.where(f"NOT {_truthy('simulated')}")

    @property
    def sql_where(self) -> str:
        return " WHERE " + " AND ".join(self.clauses) if self.clauses else ""

    def select(self, columns: list[str], order: str = "_row", limit: int | None = None) -> pd.DataFrame:
        cols = ", ".join(c for c in columns if c in self.columns)
        sql = f"SELECT {cols} FROM {self.name}{self.sql_where} ORDER BY {order}"
        if limit is not None:
            sql += f" LIMIT {int(limit)}"
        return _with_booleans(self.cur.execute(sql, self.params).df())

    def rows(self, sql: str) -> list[dict]:
        """Run `sql`, in which `{table}` stands for this filtered table, as a list of dicts."""
        sql = sql.format(table=f"(SELECT * FROM {self.name}{self.sql_where})")
        result = self.cur.execute(sql, self.params)
        names = [d[0] for d in result.description]
        return [dict(zip(names, row)) for row in result.fetchall()]

    def count_by(self, column: str) -> dict:
        """Row count per non-null value of `column`."""
        if column not in self.columns:
            return {}
        rows = self.rows(f"SELECT CAST({column} AS VARCHAR) AS v, count(*) AS n FROM {{table}} "
                         f"WHERE {column} IS NOT NULL GROUP BY v ORDER BY n DESC, v")
        return {r["v"]: int(r["n"]) for r in rows}

    def split(self) -> dict:
        """Row counts split into real and simulated."""
        if not self.n_rows or "simulated" not in self.columns:
            return {"real": 0, "simulated": 0}
        row = self.rows(f"SELECT count(*) FILTER (WHERE NOT {_truthy('simulated')}) AS real, "
                        f"count(*) FILTER (WHERE {_truthy('simulated')}) AS simulated FROM {{table}}")[0]
        return {"real": int(row["real"]), "simulated": int(row["simulated"])}

    def distinct(self, column: str) -> list[str]:
        if column not in self.columns:
            return []
        rows = self.rows(f"SELECT DISTINCT CAST({column} AS VARCHAR) AS v FROM {{table}} WHERE {column} IS NOT NULL")
        return sorted(r["v"] for r in rows)


def artifacts(sources: Sources = DEFAULT_SOURCES) -> dict:
    """Which pipeline outputs exist, and the command that produces each missing one."""
    out = {
        "evidence_table": {"present": sources.evidence_table.exists(),
                           "path": _repo_relative(sources.evidence_table),
                           "produced_by": "aree harmonize --study STUDY_ID"},
        "candidates_table": {"present": sources.candidates_table.exists(),
                             "path": _repo_relative(sources.candidates_table),
                             "produced_by": "aree build-evidence-cards"},
    }
    if sources.data_dir is not None:
        for entry in out.values():
            entry["produced_by"] = "aree export-site-data (or aree fetch-site-data for the published tables)"
    return out


# --------------------------------------------------------------------------- studies

def _study_summary(study: dict) -> dict:
    comparisons = study.get("comparisons") or []
    accessions = study.get("accessions") or {}
    availability = study.get("data_availability") or {}
    resolved = resolve_species(study.get("species")) if isinstance(study.get("species"), str) else None

    def collect(key):
        return sorted({str(c.get(key)) for c in comparisons if c.get(key)})

    return {
        "study_id": study["study_id"],
        "simulated": bool(study.get("simulated", False)),
        "species": study.get("species"),
        "species_accepted": resolved.scientific_name if resolved else study.get("species"),
        "assay_type": "|".join(study.get("assay_type") or []),
        "analysis_mode": study.get("analysis_mode"),
        "tissues": "|".join(collect("tissue")),
        "life_stages": "|".join(collect("life_stage")),
        "stressors": "|".join(collect("stressor_standardized")),
        "phenotypes": "|".join(collect("phenotype")),
        "resilience_classifications": "|".join(collect("resilience_classification")),
        "contrast_types": "|".join(collect("contrast_type")),
        "n_comparisons": len(comparisons),
        "data_status": availability.get("status"),
        "raw_data_available": availability.get("raw_data_available"),
        "processed_data_available": availability.get("processed_data_available"),
        "qc_status": study.get("qc_status"),
        "analysis_status": study.get("analysis_status"),
        "doi": accessions.get("doi"),
        "bioproject": accessions.get("bioproject"),
        "citation": study.get("citation"),
    }


def _load_studies(sources: Sources) -> list[dict]:
    return [load_yaml(p) for p in sources.study_paths()]


def study_table(sources: Sources = DEFAULT_SOURCES) -> pd.DataFrame:
    """One row per registered study, read from the study YAMLs (the source of truth)."""
    rows = [_study_summary(s) for s in _load_studies(sources)]
    return pd.DataFrame(rows, columns=STUDY_COLUMNS)


def list_studies(species: str | None = None, stressor: str | None = None,
                 phenotype: str | None = None, assay_type: str | None = None,
                 include_simulated: bool = True,
                 sources: Sources = DEFAULT_SOURCES) -> dict:
    df = study_table(sources)
    if species:
        df = df[_contains(df["species"], species) | _contains(df["species_accepted"], species)]
    for column, value in (("stressors", stressor), ("phenotypes", phenotype), ("assay_type", assay_type)):
        if value:
            df = df[_contains(df[column], value)]
    df = _filter_simulated(df, include_simulated)
    return {"n_studies": len(df), "studies": _records(df, STUDY_COLUMNS), "note": INTERPRETATION_NOTE}


def get_study(study_id: str, sources: Sources = DEFAULT_SOURCES) -> dict:
    """The full registration record for one study, plus how much evidence it contributed."""
    path = next((p for p in sources.study_paths() if p.stem == study_id), None)
    if path is None:
        known = sorted(p.stem for p in sources.study_paths())
        return {"error": f"No registered study {study_id!r}.", "known_study_ids": known}
    study = load_yaml(path)
    evidence = _Table(sources, "evidence")
    n_records = None
    if evidence.n_rows:
        n_records = evidence.where("CAST(study_id AS VARCHAR) = ?", study_id).rows(
            "SELECT count(*) AS n FROM {table}")[0]["n"]
    return {
        "study": study,
        "registration_file": _repo_relative(path),
        "n_evidence_records": n_records,
        "artifacts": artifacts(sources),
        "note": INTERPRETATION_NOTE,
    }


# --------------------------------------------------------------------------- evidence

def search_features(query: str, limit: int = 25, include_simulated: bool = True,
                    sources: Sources = DEFAULT_SOURCES) -> dict:
    """Features whose standardized/original identifier or orthogroup contains `query`,
    summarized one row per standardized feature."""
    evidence = _Table(sources, "evidence")
    if not evidence.n_rows:
        return {"n_features": 0, "features": [], "artifacts": artifacts(sources)}
    evidence.feature(query.strip(), exact=False).real_only(include_simulated)
    n = _limit(limit)
    rows = evidence.rows(f"""
        WITH hits AS (SELECT * FROM {{table}}),
        directions AS (
            SELECT feature_id_standardized, feature_type,
                   list(struct_pack(d := d, n := n) ORDER BY d) AS directions
            FROM (SELECT feature_id_standardized, feature_type,
                         CAST(molecular_direction AS VARCHAR) AS d, count(*) AS n
                  FROM hits WHERE molecular_direction IS NOT NULL GROUP BY ALL)
            GROUP BY ALL),
        features AS (
            SELECT feature_id_standardized, feature_type,
                   count(*) AS n_records,
                   count(*) FILTER (WHERE {_number('adjusted_p_value')} <= {SIGNIFICANCE_THRESHOLD})
                       AS n_significant_records,
                   count(DISTINCT study_id) AS n_studies,
                   list(DISTINCT CAST(study_id AS VARCHAR)) AS studies,
                   list(DISTINCT CAST(phenotype AS VARCHAR)) FILTER (WHERE phenotype IS NOT NULL) AS phenotypes,
                   list(DISTINCT CAST(feature_id_original AS VARCHAR))
                       FILTER (WHERE feature_id_original IS NOT NULL) AS original_ids,
                   min(CAST(orthogroup_id AS VARCHAR)) AS orthogroup_id,
                   bool_or({_truthy('simulated')}) AS any_simulated,
                   list(DISTINCT CAST(mapping_confidence AS VARCHAR))
                       FILTER (WHERE mapping_confidence IS NOT NULL) AS mapping_confidence
            FROM hits GROUP BY ALL)
        SELECT f.*, d.directions, count(*) OVER () AS n_features
        FROM features f LEFT JOIN directions d
          ON f.feature_id_standardized IS NOT DISTINCT FROM d.feature_id_standardized
         AND f.feature_type IS NOT DISTINCT FROM d.feature_type
        ORDER BY n_studies DESC, n_significant_records DESC,
                 coalesce(CAST(f.feature_id_standardized AS VARCHAR), 'None'),
                 coalesce(CAST(f.feature_type AS VARCHAR), '')
        LIMIT {n}""")
    if not rows:
        return {"n_features": 0, "features": [], "query": query}

    features = [{
        "feature_id_standardized": _clean(r["feature_id_standardized"]),
        "feature_type": _clean(r["feature_type"]),
        "original_ids": sorted(r["original_ids"] or [])[:10],
        "orthogroup_id": r["orthogroup_id"],
        "n_records": int(r["n_records"]),
        "n_significant_records": int(r["n_significant_records"]),
        "n_studies": int(r["n_studies"]),
        "studies": sorted(r["studies"] or []),
        "phenotypes": sorted(r["phenotypes"] or []),
        "directions": {d["d"]: int(d["n"]) for d in (r["directions"] or [])},
        "any_simulated": bool(r["any_simulated"]),
        "mapping_confidence": sorted(r["mapping_confidence"] or []),
    } for r in rows]
    n_features = int(rows[0]["n_features"])
    return {"query": query, "n_features": n_features, "features": features,
            "truncated": n_features > n, "note": INTERPRETATION_NOTE}


def get_evidence(feature: str | None = None, study_id: str | None = None,
                 phenotype: str | None = None, stressor: str | None = None,
                 species: str | None = None, feature_type: str | None = None,
                 tissue: str | None = None, max_adjusted_p: float | None = None,
                 include_simulated: bool = True, exact_feature: bool = True,
                 limit: int = 50, sources: Sources = DEFAULT_SOURCES) -> dict:
    """Harmonized evidence records matching every filter given, most significant first."""
    evidence = _Table(sources, "evidence")
    if not evidence.n_rows:
        return {"n_records": 0, "records": [], "artifacts": artifacts(sources)}
    if feature:
        evidence.feature(feature.strip(), exact=exact_feature)
    for column, value in (("study_id", study_id), ("phenotype", phenotype), ("stressor", stressor),
                          ("feature_type", feature_type), ("tissue", tissue)):
        if value:
            evidence.equals(column, value)
    if species:
        evidence.contains("species", species)
    if max_adjusted_p is not None:
        evidence.where(f"{_number('adjusted_p_value')} <= ?", float(max_adjusted_p))
    evidence.real_only(include_simulated)

    n = _limit(limit)
    totals = evidence.rows("SELECT count(*) AS n, count(DISTINCT study_id) AS k FROM {table}")[0]
    directions = evidence.count_by("molecular_direction")
    records = evidence.select(EVIDENCE_COLUMNS, order=f"{_number('adjusted_p_value')} ASC NULLS LAST, _row",
                              limit=n)
    return {
        "n_records": int(totals["n"]),
        "n_studies": int(totals["k"]),
        "direction_counts": directions,
        "conflicting_directions": _conflicting(directions),
        "records": _records(records, EVIDENCE_COLUMNS),
        "truncated": int(totals["n"]) > n,
        "note": INTERPRETATION_NOTE,
    }


REPLICATION_NOTE = (
    "Studies are counted by study_id. Two registered studies from one experiment (for "
    "example, sharing a control group) are not fully independent replicates; get_study "
    "shows each study's design."
)


def find_replicated_features(phenotype: str | None = None, stressor: str | None = None,
                             feature_type: str | None = None, species: str | None = None,
                             tissue: str | None = None, min_studies: int = 2,
                             max_adjusted_p: float = SIGNIFICANCE_THRESHOLD,
                             include_simulated: bool = True, limit: int = 25,
                             sources: Sources = DEFAULT_SOURCES) -> dict:
    """Features significant (adjusted p <= max_adjusted_p) in at least `min_studies`
    studies, with how their significant records agree in direction.

    Only standardized identifiers are compared: an unresolved original identifier
    cannot be matched across studies. Most studies first, then consistent before
    conflicting, then the smallest adjusted p.
    """
    evidence = _Table(sources, "evidence")
    if not evidence.n_rows:
        return {"n_features": 0, "features": [], "artifacts": artifacts(sources)}
    for column, value in (("phenotype", phenotype), ("stressor", stressor),
                          ("feature_type", feature_type), ("tissue", tissue)):
        if value:
            evidence.equals(column, value)
    if species:
        evidence.contains("species", species)
    evidence.where(f"{_number('adjusted_p_value')} <= ?", float(max_adjusted_p))
    evidence.where("feature_id_standardized IS NOT NULL")
    evidence.real_only(include_simulated)

    n = _limit(limit)
    up = "('up', 'hyper')"
    down = "('down', 'hypo')"
    direction = "lower(CAST(molecular_direction AS VARCHAR))"
    rows = evidence.rows(f"""
        WITH hits AS (SELECT * FROM {{table}}),
        directions AS (
            SELECT feature_id_standardized, feature_type,
                   list(struct_pack(d := d, n := n) ORDER BY d) AS directions
            FROM (SELECT feature_id_standardized, feature_type,
                         CAST(molecular_direction AS VARCHAR) AS d, count(*) AS n
                  FROM hits WHERE molecular_direction IS NOT NULL GROUP BY ALL)
            GROUP BY ALL),
        features AS (
            SELECT feature_id_standardized, feature_type,
                   count(DISTINCT study_id) AS n_studies,
                   count(*) AS n_significant_records,
                   count(DISTINCT study_id) FILTER (WHERE {direction} IN {up}) AS n_studies_up,
                   count(DISTINCT study_id) FILTER (WHERE {direction} IN {down}) AS n_studies_down,
                   min({_number('adjusted_p_value')}) AS min_adjusted_p,
                   list(DISTINCT CAST(study_id AS VARCHAR)) AS studies,
                   list(DISTINCT CAST(phenotype AS VARCHAR)) FILTER (WHERE phenotype IS NOT NULL) AS phenotypes,
                   list(DISTINCT CAST(stressor AS VARCHAR)) FILTER (WHERE stressor IS NOT NULL) AS stressors,
                   list(DISTINCT CAST(tissue AS VARCHAR)) FILTER (WHERE tissue IS NOT NULL) AS tissues,
                   list(DISTINCT CAST(mapping_confidence AS VARCHAR))
                       FILTER (WHERE mapping_confidence IS NOT NULL) AS mapping_confidence,
                   min(CAST(orthogroup_id AS VARCHAR)) AS orthogroup_id,
                   bool_or({_truthy('simulated')}) AS any_simulated
            FROM hits GROUP BY ALL
            HAVING count(DISTINCT study_id) >= {int(min_studies)})
        SELECT f.*, d.directions,
               (f.n_studies_up > 0 AND f.n_studies_down > 0) AS conflicting_directions,
               count(*) OVER () AS n_features,
               count(*) FILTER (WHERE f.n_studies_up > 0 AND f.n_studies_down > 0) OVER () AS n_conflicting
        FROM features f LEFT JOIN directions d
          ON f.feature_id_standardized = d.feature_id_standardized
         AND f.feature_type IS NOT DISTINCT FROM d.feature_type
        ORDER BY n_studies DESC, conflicting_directions, min_adjusted_p,
                 CAST(f.feature_id_standardized AS VARCHAR), coalesce(CAST(f.feature_type AS VARCHAR), '')
        LIMIT {n}""")
    criteria = {"min_studies": int(min_studies), "max_adjusted_p": float(max_adjusted_p)}
    if not rows:
        return {"n_features": 0, "n_conflicting": 0, "features": [], "criteria": criteria,
                "note": f"{REPLICATION_NOTE} {INTERPRETATION_NOTE}"}

    features = [{
        "feature_id_standardized": r["feature_id_standardized"],
        "feature_type": r["feature_type"],
        "orthogroup_id": r["orthogroup_id"],
        "n_studies": int(r["n_studies"]),
        "n_significant_records": int(r["n_significant_records"]),
        "n_studies_up": int(r["n_studies_up"]),
        "n_studies_down": int(r["n_studies_down"]),
        "conflicting_directions": bool(r["conflicting_directions"]),
        "directions": {d["d"]: int(d["n"]) for d in (r["directions"] or [])},
        "min_adjusted_p": _clean(r["min_adjusted_p"]),
        "studies": sorted(r["studies"] or []),
        "phenotypes": sorted(r["phenotypes"] or []),
        "stressors": sorted(r["stressors"] or []),
        "tissues": sorted(r["tissues"] or []),
        "mapping_confidence": sorted(r["mapping_confidence"] or []),
        "any_simulated": bool(r["any_simulated"]),
    } for r in rows]
    n_features = int(rows[0]["n_features"])
    return {"n_features": n_features, "n_conflicting": int(rows[0]["n_conflicting"]),
            "features": features, "truncated": n_features > n, "criteria": criteria,
            "note": f"{REPLICATION_NOTE} {INTERPRETATION_NOTE}"}


# --------------------------------------------------------------------------- candidates

def list_candidates(phenotype: str | None = None, tier: str | None = None,
                    evidence_class: str | None = None, feature: str | None = None,
                    species_taxid: int | None = None, feature_type: str | None = None,
                    include_simulated: bool = True, limit: int = 25,
                    contrast_type: str | None = None,
                    sources: Sources = DEFAULT_SOURCES) -> dict:
    """Ranked candidate biomarkers from `aree build-evidence-cards`, best score first within tier order.

    contrast_type="baseline_phenotype_difference" returns candidates resting on a
    difference between resilient and susceptible groups under unstressed conditions.
    """
    candidates = _Table(sources, "candidates")
    if not candidates.n_rows:
        return {"n_candidates": 0, "candidates": [], "artifacts": artifacts(sources)}
    for column, value in (("phenotype", phenotype), ("tier", tier), ("evidence_class", evidence_class),
                          ("feature_type", feature_type), ("contrast_type", contrast_type)):
        if value:
            candidates.equals(column, value)
    if feature:
        candidates.contains("feature_id_standardized", feature.strip())
    if species_taxid is not None:
        candidates.where(f"{_number('species_taxid')} = ?", int(species_taxid))
    candidates.real_only(include_simulated)

    n = _limit(limit)
    tiers = candidates.count_by("tier")
    n_candidates = candidates.rows("SELECT count(*) AS n FROM {table}")[0]["n"]
    records = _records(candidates.select(CANDIDATE_COLUMNS, limit=n), CANDIDATE_COLUMNS)
    for rec in records:
        if rec.get("card_file"):
            rec["card_file"] = _repo_relative(rec["card_file"])
    return {"n_candidates": int(n_candidates), "tier_counts": tiers,
            "candidates": records, "truncated": n_candidates > n, "note": INTERPRETATION_NOTE}


def _card_markdown(card_file, sources: Sources) -> tuple[str | None, str | None]:
    """(where the card was found, its markdown), or (None, None) when no card was rendered."""
    if not isinstance(card_file, str) or not card_file:
        return None, None
    name = Path(card_file).name
    if sources.cards_table is not None:
        if not sources.cards_table.exists():
            return None, None
        row = sources.db().execute(
            f"SELECT markdown FROM read_parquet({_sql_literal(sources.cards_table)}) WHERE card_file = ?",
            [name]).fetchone()
        return (name, row[0]) if row else (None, None)
    path = Path(card_file)
    if not path.exists():  # written on another machine: resolve by name here
        path = sources.cards_dir / name
    return (_repo_relative(path), path.read_text()) if path.exists() else (None, None)


def get_evidence_card(feature: str, phenotype: str | None = None,
                      include_simulated: bool = True,
                      sources: Sources = DEFAULT_SOURCES) -> dict:
    """The rendered markdown evidence card(s) for a candidate feature.

    Cards are only written for candidates with a significant signal; for any
    other candidate the ranked row is returned with an explanation instead.
    """
    candidates = _Table(sources, "candidates")
    if not candidates.n_rows:
        return {"n_cards": 0, "cards": [], "artifacts": artifacts(sources)}
    candidates.where(f"{_text('feature_id_standardized')} = lower(?)", feature.strip())
    if phenotype:
        candidates.equals("phenotype", phenotype)
    candidates.real_only(include_simulated)
    rows = candidates.select(CANDIDATE_COLUMNS)
    if rows.empty:
        return {"n_cards": 0, "cards": [],
                "message": f"No ranked candidate {feature!r}"
                           + (f" for phenotype {phenotype!r}" if phenotype else "")
                           + ". Try search_features to find the standardized identifier."}

    cards = []
    for row in rows.to_dict("records"):
        entry = {"candidate": {k: _clean(row.get(k)) for k in CANDIDATE_COLUMNS if k in row}}
        found, markdown = _card_markdown(row.get("card_file"), sources)
        if markdown is not None:
            entry["card_file"] = found
            entry["markdown"] = markdown
        else:
            entry["markdown"] = None
            entry["message"] = ("No card was rendered for this candidate: cards are written only for "
                                "candidates with a significant signal (aree build-evidence-cards --all-cards "
                                "renders every one).")
        cards.append(entry)
    return {"n_cards": len(cards), "cards": cards, "note": INTERPRETATION_NOTE}


# --------------------------------------------------------------------------- vocabularies + overview

VOCABULARIES = {
    "phenotype": "phenotype_ontology",
    "stressor": "stressor_ontology",
    "assay_type": "assay_types",
    "feature_type": "feature_types",
    "tissue": "tissue_types",
    "life_stage": "life_stages",
    "mapping_confidence": "mapping_confidence",
    "quality_flag": "quality_flags",
    "species": "species",
}


def list_vocabulary(name: str) -> dict:
    """Term ids and labels for a controlled vocabulary, so free text can be mapped to filter values."""
    stem = VOCABULARIES.get(name, name)
    try:
        vocab = load_vocab(stem)
    except FileNotFoundError:
        return {"error": f"Unknown vocabulary {name!r}.", "available": sorted(VOCABULARIES)}
    terms = []
    for tid, term in vocab.items():
        terms.append({
            "id": tid,
            "label": term.get("label") or term.get("scientific_name"),
            "definition": term.get("definition") or term.get("description"),
        })
    return {"vocabulary": stem, "n_terms": len(terms), "terms": terms}


def overview(sources: Sources = DEFAULT_SOURCES) -> dict:
    """Counts of what AREE currently holds, split real vs simulated."""
    studies = study_table(sources)
    evidence = _Table(sources, "evidence")
    candidates = _Table(sources, "candidates")

    if studies.empty:
        study_split = {"real": 0, "simulated": 0}
    else:
        sim = _as_bool_series(studies["simulated"])
        study_split = {"real": int((~sim).sum()), "simulated": int(sim.sum())}

    out = {
        "studies": study_split,
        "evidence_records": evidence.split(),
        "candidates": candidates.split(),
        "artifacts": artifacts(sources),
        "note": INTERPRETATION_NOTE,
    }
    if evidence.n_rows:
        out["phenotypes_with_evidence"] = evidence.distinct("phenotype")
        out["stressors_with_evidence"] = evidence.distinct("stressor")
        out["feature_types"] = evidence.distinct("feature_type")
        out["species_with_evidence"] = evidence.distinct("species")
    if candidates.n_rows and "tier" in candidates.columns:
        out["candidate_tiers"] = candidates.count_by("tier")
    return out
