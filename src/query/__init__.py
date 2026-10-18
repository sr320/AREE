"""Read-only query layer over AREE's registry and pipeline outputs.

One set of functions answers "what does AREE hold?" questions — studies,
harmonized evidence, ranked candidates, evidence cards, vocabularies — for the
MCP server (`query.mcp_server`) and the static-site export (`query.site_export`).

It reads what the `aree` pipeline has already written; it never recomputes a
meta-analysis or writes anything. Every function returns plain JSON-serializable
dicts/lists, and every record keeps its provenance (study_id, evidence_id,
source_file, workflow_version, date_generated) and its `simulated` flag, so a
caller — human or language model — can always say where a number came from.

Missing pipeline outputs are reported as missing (`artifacts`), never silently
treated as "no evidence".
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from common import REPO_ROOT, REPORTS_DIR, STUDIES_DIR, load_vocab, load_yaml, resolve_species

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


@dataclass
class Sources:
    """Where the query layer reads from; tests point these at a temp directory."""

    reports_dir: Path = REPORTS_DIR
    studies_dir: Path = STUDIES_DIR
    _cache: dict = field(default_factory=dict, repr=False)

    @property
    def evidence_table(self) -> Path:
        return self.reports_dir / "evidence" / "evidence_table.tsv"

    @property
    def candidates_table(self) -> Path:
        return self.reports_dir / "evidence_cards" / "candidates.tsv"

    @property
    def cards_dir(self) -> Path:
        return self.reports_dir / "evidence_cards"

    def _read_tsv(self, path: Path) -> pd.DataFrame:
        if not path.exists():
            return pd.DataFrame()
        key = (str(path), path.stat().st_mtime_ns)
        if key not in self._cache:
            df = pd.read_csv(path, sep="\t", dtype=str)
            for col in NUMERIC_COLUMNS & set(df.columns):
                df[col] = pd.to_numeric(df[col], errors="coerce")
            self._cache[key] = df
        return self._cache[key]

    def evidence(self) -> pd.DataFrame:
        return self._read_tsv(self.evidence_table)

    def candidates(self) -> pd.DataFrame:
        return self._read_tsv(self.candidates_table)


DEFAULT_SOURCES = Sources()


# --------------------------------------------------------------------------- helpers

def _clean(value):
    """NaN/inf -> None, numpy scalars -> Python, so results serialize as JSON."""
    if value is None:
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


def _contains(series: pd.Series, needle: str) -> pd.Series:
    return series.astype(str).str.lower().str.contains(needle.lower(), na=False, regex=False)


def _filter_simulated(df: pd.DataFrame, include_simulated: bool) -> pd.DataFrame:
    if include_simulated or "simulated" not in df.columns or df.empty:
        return df
    return df[~_as_bool_series(df["simulated"])]


def _repo_relative(path) -> str:
    try:
        return str(Path(path).resolve().relative_to(REPO_ROOT))
    except (ValueError, OSError):
        return str(path)


def artifacts(sources: Sources = DEFAULT_SOURCES) -> dict:
    """Which pipeline outputs exist, and the command that produces each missing one."""
    return {
        "evidence_table": {"present": sources.evidence_table.exists(),
                           "path": _repo_relative(sources.evidence_table),
                           "produced_by": "aree harmonize --study STUDY_ID"},
        "candidates_table": {"present": sources.candidates_table.exists(),
                             "path": _repo_relative(sources.candidates_table),
                             "produced_by": "aree build-evidence-cards"},
    }


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
    paths = sorted(p for p in sources.studies_dir.glob("*.yaml") if not p.name.startswith("_"))
    return [load_yaml(p) for p in paths]


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
    path = sources.studies_dir / f"{study_id}.yaml"
    if not path.exists():
        known = sorted(p.stem for p in sources.studies_dir.glob("*.yaml") if not p.name.startswith("_"))
        return {"error": f"No registered study {study_id!r}.", "known_study_ids": known}
    study = load_yaml(path)
    evidence = sources.evidence()
    n_records = int((evidence["study_id"] == study_id).sum()) if not evidence.empty else None
    return {
        "study": study,
        "registration_file": _repo_relative(path),
        "n_evidence_records": n_records,
        "artifacts": artifacts(sources),
        "note": INTERPRETATION_NOTE,
    }


# --------------------------------------------------------------------------- evidence

def _feature_mask(df: pd.DataFrame, feature: str, exact: bool) -> pd.Series:
    cols = [c for c in ("feature_id_standardized", "feature_id_original", "orthogroup_id") if c in df.columns]
    mask = pd.Series(False, index=df.index)
    for col in cols:
        if exact:
            mask |= df[col].astype(str).str.lower() == feature.lower()
        else:
            mask |= _contains(df[col], feature)
    return mask


def search_features(query: str, limit: int = 25, include_simulated: bool = True,
                    sources: Sources = DEFAULT_SOURCES) -> dict:
    """Features whose standardized/original identifier or orthogroup contains `query`,
    summarized one row per standardized feature."""
    evidence = sources.evidence()
    if evidence.empty:
        return {"n_features": 0, "features": [], "artifacts": artifacts(sources)}
    hits = _filter_simulated(evidence[_feature_mask(evidence, query.strip(), exact=False)], include_simulated)
    if hits.empty:
        return {"n_features": 0, "features": [], "query": query}

    adj = pd.to_numeric(hits["adjusted_p_value"], errors="coerce")
    hits = hits.assign(_sig=adj <= 0.05)
    rows = []
    for (feature, ftype), grp in hits.groupby(["feature_id_standardized", "feature_type"], dropna=False):
        directions = grp["molecular_direction"].astype(str).value_counts().to_dict()
        rows.append({
            "feature_id_standardized": _clean(feature),
            "feature_type": _clean(ftype),
            "original_ids": sorted(grp["feature_id_original"].astype(str).unique())[:10],
            "orthogroup_id": _clean(grp["orthogroup_id"].dropna().iloc[0]) if grp["orthogroup_id"].notna().any() else None,
            "n_records": len(grp),
            "n_significant_records": int(grp["_sig"].sum()),
            "n_studies": int(grp["study_id"].nunique()),
            "studies": sorted(grp["study_id"].astype(str).unique()),
            "phenotypes": sorted(grp["phenotype"].dropna().astype(str).unique()),
            "directions": {k: int(v) for k, v in directions.items()},
            "any_simulated": bool(_as_bool_series(grp["simulated"]).any()),
            "mapping_confidence": sorted(grp["mapping_confidence"].dropna().astype(str).unique()),
        })
    rows.sort(key=lambda r: (-r["n_studies"], -r["n_significant_records"], str(r["feature_id_standardized"])))
    n = _limit(limit)
    return {"query": query, "n_features": len(rows), "features": rows[:n],
            "truncated": len(rows) > n, "note": INTERPRETATION_NOTE}


def get_evidence(feature: str | None = None, study_id: str | None = None,
                 phenotype: str | None = None, stressor: str | None = None,
                 species: str | None = None, feature_type: str | None = None,
                 tissue: str | None = None, max_adjusted_p: float | None = None,
                 include_simulated: bool = True, exact_feature: bool = True,
                 limit: int = 50, sources: Sources = DEFAULT_SOURCES) -> dict:
    """Harmonized evidence records matching every filter given, most significant first."""
    df = sources.evidence()
    if df.empty:
        return {"n_records": 0, "records": [], "artifacts": artifacts(sources)}
    if feature:
        df = df[_feature_mask(df, feature.strip(), exact=exact_feature)]
    for column, value in (("study_id", study_id), ("phenotype", phenotype), ("stressor", stressor),
                          ("feature_type", feature_type), ("tissue", tissue)):
        if value:
            df = df[df[column].astype(str).str.lower() == value.lower()]
    if species:
        df = df[_contains(df["species"], species)]
    if max_adjusted_p is not None:
        df = df[pd.to_numeric(df["adjusted_p_value"], errors="coerce") <= float(max_adjusted_p)]
    df = _filter_simulated(df, include_simulated)
    df = df.assign(_p=pd.to_numeric(df["adjusted_p_value"], errors="coerce")).sort_values("_p", na_position="last")

    n = _limit(limit)
    directions = df["molecular_direction"].astype(str).value_counts().to_dict() if len(df) else {}
    return {
        "n_records": len(df),
        "n_studies": int(df["study_id"].nunique()) if len(df) else 0,
        "direction_counts": {k: int(v) for k, v in directions.items()},
        "conflicting_directions": len({d for d in directions if d in ("up", "down")}) > 1,
        "records": _records(df.head(n), EVIDENCE_COLUMNS),
        "truncated": len(df) > n,
        "note": INTERPRETATION_NOTE,
    }


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
    df = sources.candidates()
    if df.empty:
        return {"n_candidates": 0, "candidates": [], "artifacts": artifacts(sources)}
    for column, value in (("phenotype", phenotype), ("tier", tier), ("evidence_class", evidence_class),
                          ("feature_type", feature_type), ("contrast_type", contrast_type)):
        if value and column in df.columns:
            df = df[df[column].astype(str).str.lower() == value.lower()]
    if feature:
        df = df[_contains(df["feature_id_standardized"], feature.strip())]
    if species_taxid is not None:
        df = df[pd.to_numeric(df["species_taxid"], errors="coerce") == int(species_taxid)]
    df = _filter_simulated(df, include_simulated)

    n = _limit(limit)
    records = _records(df.head(n), CANDIDATE_COLUMNS)
    for rec in records:
        if rec.get("card_file"):
            rec["card_file"] = _repo_relative(rec["card_file"])
    tiers = df["tier"].astype(str).value_counts().to_dict() if "tier" in df.columns else {}
    return {"n_candidates": len(df), "tier_counts": {k: int(v) for k, v in tiers.items()},
            "candidates": records, "truncated": len(df) > n, "note": INTERPRETATION_NOTE}


def get_evidence_card(feature: str, phenotype: str | None = None,
                      include_simulated: bool = True,
                      sources: Sources = DEFAULT_SOURCES) -> dict:
    """The rendered markdown evidence card(s) for a candidate feature.

    Cards are only written for candidates with a significant signal; for any
    other candidate the ranked row is returned with an explanation instead.
    """
    df = sources.candidates()
    if df.empty:
        return {"n_cards": 0, "cards": [], "artifacts": artifacts(sources)}
    df = df[df["feature_id_standardized"].astype(str).str.lower() == feature.strip().lower()]
    if phenotype:
        df = df[df["phenotype"].astype(str).str.lower() == phenotype.lower()]
    df = _filter_simulated(df, include_simulated)
    if df.empty:
        return {"n_cards": 0, "cards": [],
                "message": f"No ranked candidate {feature!r}"
                           + (f" for phenotype {phenotype!r}" if phenotype else "")
                           + ". Try search_features to find the standardized identifier."}

    cards = []
    for row in df.to_dict("records"):
        card_file = row.get("card_file")
        path = None
        if isinstance(card_file, str) and card_file:
            path = Path(card_file)
            if not path.exists():  # written on another machine: resolve by name here
                path = sources.cards_dir / path.name
        entry = {"candidate": {k: _clean(row.get(k)) for k in CANDIDATE_COLUMNS if k in row}}
        if path is not None and path.exists():
            entry["card_file"] = _repo_relative(path)
            entry["markdown"] = path.read_text()
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
    evidence = sources.evidence()
    candidates = sources.candidates()

    def split(df):
        if df.empty or "simulated" not in df.columns:
            return {"real": 0, "simulated": 0}
        sim = _as_bool_series(df["simulated"])
        return {"real": int((~sim).sum()), "simulated": int(sim.sum())}

    out = {
        "studies": split(studies),
        "evidence_records": split(evidence),
        "candidates": split(candidates),
        "artifacts": artifacts(sources),
        "note": INTERPRETATION_NOTE,
    }
    if not evidence.empty:
        out["phenotypes_with_evidence"] = sorted(evidence["phenotype"].dropna().astype(str).unique())
        out["stressors_with_evidence"] = sorted(evidence["stressor"].dropna().astype(str).unique())
        out["feature_types"] = sorted(evidence["feature_type"].dropna().astype(str).unique())
        out["species_with_evidence"] = sorted(evidence["species"].dropna().astype(str).unique())
    if not candidates.empty and "tier" in candidates.columns:
        out["candidate_tiers"] = {k: int(v) for k, v in candidates["tier"].astype(str).value_counts().items()}
    return out
