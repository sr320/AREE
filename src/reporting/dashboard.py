"""Progress-and-findings dashboard data for the GitHub Pages site.

`build_dashboard_data()` reads what the repository and the `aree` pipeline have
actually produced — registry YAMLs, the evidence table, meta-analysis tables,
ranked candidates, manifests, the test suite, the git log, and the candid
implementation-status page — and writes one JSON document that
`docs/index.qmd` (the site home page) renders.

Everything is derived; nothing is typed in by hand. Missing pipeline outputs are
reported as missing (`artifacts_present`), not silently zeroed, so a dashboard
built without running the pipeline says so on the page.

Real and simulated evidence are kept apart in every table here, exactly as they
are in pooling and ranking.
"""
from __future__ import annotations

import datetime as _dt
import json
import re
import subprocess
from collections import Counter
from pathlib import Path

import pandas as pd
import yaml

from common import (
    REPO_ROOT,
    REPORTS_DIR,
    load_vocab,
    registered_crosswalks,
    resolve_species,
    species_term_for_taxid,
    study_files,
)
from prioritize.rank import DEFAULT_CONTRAST_TYPE

from .species import species_identity_basis, species_summary

EVIDENCE_TABLE = REPORTS_DIR / "evidence" / "evidence_table.tsv"
CANDIDATES_TABLE = REPORTS_DIR / "evidence_cards" / "candidates.tsv"
CROSS_SPECIES_TABLE = REPORTS_DIR / "evidence_cards" / "cross_species" / "candidates.tsv"
META_ANALYSIS_DIR = REPORTS_DIR / "meta_analysis"
MANIFESTS_DIR = REPORTS_DIR / "manifests"
TOP_CANDIDATES_MD = REPORTS_DIR / "top_candidates_summary.md"
IMPLEMENTATION_STATUS_MD = REPO_ROOT / "docs" / "implementation_status.md"


def _built_crosswalks() -> list[Path]:
    """Every real crosswalk registered in genome_assemblies.yaml that exists here."""
    return [entry["path"] for entry in registered_crosswalks() if entry["path"].exists()]
TESTS_DIR = REPO_ROOT / "tests"
DEFAULT_OUT_PATH = REPO_ROOT / "docs" / "dashboard" / "data.json"
CARDS_SUBDIR = "cards"  # evidence cards for the top candidates, copied beside data.json
REPO_URL = "https://github.com/sr320/AREE"

SIGNIFICANCE_THRESHOLD = 0.05
TOP_N_CANDIDATES = 15
# The overall top list is dominated by whichever pool has the most studies, so
# each phenotype also gets its own short list.
TOP_N_PER_PHENOTYPE = 5
# Significant predictive baseline candidates shown on the home page: features
# that differ between resilient and susceptible groups under unstressed
# conditions (docs/contrast_types.md).
TOP_N_PREDICTIVE_BASELINE = 25
TOP_N_CROSS_SPECIES = 15
# WP6 dashboard gate: a phenotype's cross-species pooled results are shown only
# once at least this many species each have a within-species pool with k >= 2.
GATE_MIN_SPECIES_K2 = 2


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #


def _git(*args: str) -> str:
    try:
        out = subprocess.run(
            ["git", *args], cwd=REPO_ROOT, capture_output=True, text=True, check=True
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""
    return out.stdout.strip()


def _as_list(value) -> list:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _bool(value) -> bool:
    if isinstance(value, str):
        return value.strip().lower() == "true"
    return bool(value)


def _repo_relative(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:  # outputs redirected elsewhere (tests, --out)
        return path.name


# --------------------------------------------------------------------------- #
# Registry
# --------------------------------------------------------------------------- #


def _drop_simulated(df: pd.DataFrame, include_simulated: bool) -> pd.DataFrame:
    """Normalize the `simulated` column and, unless asked to keep them, drop simulated rows."""
    df["simulated"] = df["simulated"].map(_bool)
    return df if include_simulated else df[~df["simulated"]].reset_index(drop=True)


def _load_studies(include_simulated: bool = False) -> list[dict]:
    studies = []
    for path in study_files(include_demo=include_simulated):
        with open(path) as fh:
            doc = yaml.safe_load(fh)
        if isinstance(doc, dict) and doc.get("study_id"):
            doc["_path"] = str(path.relative_to(REPO_ROOT))
            studies.append(doc)
    return studies


def _primary_accession(study: dict) -> str:
    acc = study.get("accessions") or {}
    for key in ("bioproject", "geo", "sra", "ena", "proteomexchange", "doi", "other"):
        if acc.get(key):
            return f"{key}:{acc[key]}" if key != "doi" else f"doi:{acc[key]}"
    return ""


def _short_citation(citation: str | None) -> str:
    if not citation:
        return ""
    # "Author A, Author B (2026). Title..." -> "Author A, Author B (2026)"
    m = re.match(r"^(.*?\(\d{4}[a-z]?\))", citation)
    return m.group(1) if m else citation[:80]


def _study_taxid(study: dict):
    resolved = resolve_species(study.get("species"))
    return resolved.ncbi_taxid if resolved else None


def _study_canonical(study: dict):
    resolved = resolve_species(study.get("species"))
    return resolved.scientific_name if resolved else study.get("species")


def _study_complex(study: dict):
    resolved = resolve_species(study.get("species"))
    return resolved.complex if resolved else None


def _study_rows(studies: list[dict], evidence: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    records_by_comparison: Counter = Counter()
    if len(evidence):
        records_by_comparison.update(
            zip(evidence["study_id"], evidence["comparison_id"])
        )

    study_rows, comparison_rows = [], []
    for study in studies:
        sid = study["study_id"]
        simulated = _bool(study.get("simulated", False))
        comparisons = study.get("comparisons") or []
        n_harmonized = 0
        n_with_results = 0
        for comp in comparisons:
            cid = comp.get("comparison_id")
            n_records = records_by_comparison.get((sid, cid), 0)
            results_file = comp.get("results_file")
            results_present = bool(results_file) and (REPO_ROOT / results_file).exists()
            manifest_present = (MANIFESTS_DIR / f"{sid}_{cid}_manifest.json").exists()
            n_harmonized += int(n_records > 0)
            n_with_results += int(results_present)
            comparison_rows.append({
                "study_id": sid,
                "comparison_id": cid,
                "simulated": simulated,
                "tissue": comp.get("tissue"),
                "life_stage": comp.get("life_stage"),
                "stressor": comp.get("stressor_standardized"),
                "stressor_original": comp.get("stressor_original"),
                "phenotype": comp.get("phenotype"),
                "resilience_classification": comp.get("resilience_classification"),
                "sample_size": comp.get("sample_size"),
                "results_file_present": results_present,
                "harmonized": n_records > 0,
                "n_evidence_records": n_records,
                "manifest_present": manifest_present,
                "meta_analysis_primary": bool(comp.get("meta_analysis_primary", False)),
                "n_excluded_samples": len(comp.get("excluded_samples") or []),
            })

        availability = study.get("data_availability") or {}
        study_rows.append({
            "study_id": sid,
            "simulated": simulated,
            "species": study.get("species"),
            "species_taxid": _study_taxid(study),
            # The vocabulary's accepted name, so a synonym (Crassostrea gigas)
            # groups with its species (Magallana gigas) on the page.
            "species_canonical": _study_canonical(study),
            "species_complex": _study_complex(study),
            "species_identity_basis": species_identity_basis(study),
            "assay_types": _as_list(study.get("assay_type")),
            "analysis_mode": study.get("analysis_mode"),
            "analysis_status": study.get("analysis_status"),
            "qc_status": study.get("qc_status"),
            "data_status": availability.get("status"),
            "raw_data_available": availability.get("raw_data_available"),
            "processed_data_available": availability.get("processed_data_available"),
            "genome_assembly": study.get("genome_assembly"),
            "accession": _primary_accession(study),
            "doi": (study.get("accessions") or {}).get("doi"),
            "citation": _short_citation(study.get("citation")),
            "stressors": sorted({c.get("stressor_standardized") for c in comparisons
                                 if c.get("stressor_standardized")}),
            "phenotypes": sorted({c.get("phenotype") for c in comparisons if c.get("phenotype")}),
            "resilience_classifications": sorted({c.get("resilience_classification")
                                                  for c in comparisons
                                                  if c.get("resilience_classification")}),
            "tissues": sorted({c.get("tissue") for c in comparisons if c.get("tissue")}),
            "life_stages": sorted({c.get("life_stage") for c in comparisons if c.get("life_stage")}),
            "n_comparisons": len(comparisons),
            "n_comparisons_with_results": n_with_results,
            "n_comparisons_harmonized": n_harmonized,
            "n_evidence_records": int(sum(
                records_by_comparison.get((sid, c.get("comparison_id")), 0) for c in comparisons
            )),
            "quality_flags": _as_list(study.get("quality_flags")),
            "date_registered": str((study.get("provenance") or {}).get("date_registered") or ""),
            "registry_file": study["_path"],
        })
    return study_rows, comparison_rows


# --------------------------------------------------------------------------- #
# Evidence, meta-analysis, candidates
# --------------------------------------------------------------------------- #


def _load_evidence(include_simulated: bool = False) -> pd.DataFrame:
    if not EVIDENCE_TABLE.exists():
        return pd.DataFrame()
    df = pd.read_csv(EVIDENCE_TABLE, sep="\t", low_memory=False)
    return _drop_simulated(df, include_simulated)


def _evidence_summary(evidence: pd.DataFrame) -> dict:
    if not len(evidence):
        return {"by_study": [], "mapping_confidence": [], "by_feature_type": []}
    adj = pd.to_numeric(evidence["adjusted_p_value"], errors="coerce")
    sig = adj <= SIGNIFICANCE_THRESHOLD
    unresolved = evidence["mapping_confidence"].eq("unresolved")
    by_study = []
    for (sid, simulated), grp in evidence.groupby(["study_id", "simulated"], sort=True):
        idx = grp.index
        by_study.append({
            "study_id": sid,
            "simulated": bool(simulated),
            "n_records": int(len(grp)),
            "n_significant": int(sig.loc[idx].sum()),
            "n_unresolved": int(unresolved.loc[idx].sum()),
            "n_poolable": int(pd.to_numeric(grp["standard_error"], errors="coerce").notna().sum()),
            "feature_types": sorted(grp["feature_type"].dropna().unique().tolist()),
        })
    mapping = [
        {"simulated": bool(s), "mapping_confidence": mc, "n": int(n)}
        for (s, mc), n in evidence.groupby(["simulated", "mapping_confidence"]).size().items()
    ]
    by_feature_type = [
        {"simulated": bool(s), "feature_type": ft, "n": int(n)}
        for (s, ft), n in evidence.groupby(["simulated", "feature_type"]).size().items()
    ]
    return {"by_study": by_study, "mapping_confidence": mapping, "by_feature_type": by_feature_type}


def gene_annotation_frame() -> pd.DataFrame:
    """ncbi_gene_id, gene_symbol, gene_description for every gene in the built
    crosswalks, one row per GeneID. Empty when no crosswalk is built (e.g. CI
    without the reference download)."""
    frames = [
        pd.read_csv(cw, sep="\t", comment="#", dtype=str,
                    usecols=["ncbi_gene_id", "gene_symbol", "gene_description"])
        for cw in _built_crosswalks()
    ]
    if not frames:
        return pd.DataFrame(columns=["ncbi_gene_id", "gene_symbol", "gene_description"])
    return pd.concat(frames, ignore_index=True).dropna(subset=["ncbi_gene_id"]).drop_duplicates("ncbi_gene_id")


def _gene_annotations(gene_ids: list[str]) -> dict[str, dict]:
    """Symbol and description for real NCBI GeneIDs, from the committed crosswalks.

    NCBI GeneIDs are unique across all taxa, so the crosswalks of every
    registered species can be searched together without ambiguity.
    """
    crosswalks = _built_crosswalks()
    if not crosswalks or not gene_ids:
        return {}
    wanted = set(map(str, gene_ids))
    out: dict[str, dict] = {}
    for crosswalk in crosswalks:
        for chunk in pd.read_csv(crosswalk, sep="\t", comment="#", dtype=str,
                                 chunksize=20000):
            hit = chunk[chunk["ncbi_gene_id"].isin(wanted)]
            for _, row in hit.iterrows():
                out[row["ncbi_gene_id"]] = {
                    "gene_symbol": row.get("gene_symbol") or "",
                    "gene_description": row.get("gene_description") or "",
                    "gene_type": row.get("gene_type") or "",
                }
            if len(out) >= len(wanted):
                return out
    return out


def _candidate_summary(phenotype_labels: dict, cards_dir: Path | None = None,
                       include_simulated: bool = False) -> dict:
    empty = {"by_tier": [], "by_evidence_class": [], "top_real": [], "top_simulated": [],
             "top_real_by_phenotype": {}, "top_simulated_by_phenotype": {}, "n_total": 0,
             "n_real": 0, "n_simulated": 0, "predictive_baseline": _empty_predictive_baseline()}
    if not CANDIDATES_TABLE.exists():
        return empty
    df = pd.read_csv(CANDIDATES_TABLE, sep="\t", low_memory=False)
    df = _drop_simulated(df, include_simulated)
    if not len(df):
        return empty
    by_tier = [
        {"simulated": bool(s), "tier": t, "n": int(n)}
        for (s, t), n in df.groupby(["simulated", "tier"]).size().items()
    ]
    by_class = [
        {"simulated": bool(s), "evidence_class": c, "n": int(n)}
        for (s, c), n in df.groupby(["simulated", "evidence_class"]).size().items()
    ]

    from prioritize.rank import TIER_ORDER  # local import keeps module import cheap

    tier_rank = {t: i for i, t in enumerate(TIER_ORDER)}
    df["_tier_rank"] = df["tier"].map(tier_rank).fillna(len(TIER_ORDER))
    ranked = df.sort_values(["_tier_rank", "score"], ascending=[True, False])

    cols = ["feature_id_standardized", "feature_type", "phenotype", "tier", "evidence_class",
            "context_replication", "contrast_type", "is_predictive_baseline",
            "score", "pooling_status", "k_studies", "studies", "total_sample_size",
            "pooled_effect", "adjusted_p_value", "direction_consistency", "i_squared",
            "n_supporting_layers", "mapping_confidences", "simulated", "species_taxid",
            "orthogroup_id", "card_file"]
    cols = [c for c in cols if c in ranked.columns]

    def _top(frame: pd.DataFrame, n: int = TOP_N_CANDIDATES) -> list[dict]:
        rows = frame.head(n)[cols].to_dict(orient="records")
        for r in rows:
            r["feature_id_standardized"] = str(r["feature_id_standardized"])
            if "species_taxid" in r:
                r["species_taxid"] = _taxid(r["species_taxid"])
                r["species"] = _species_label(r["species_taxid"])
            r["phenotype_label"] = phenotype_labels.get(r.get("phenotype"), r.get("phenotype"))
            if "is_predictive_baseline" in r:
                r["is_predictive_baseline"] = str(r["is_predictive_baseline"]).strip().lower() == "true"
            for k, v in list(r.items()):
                if isinstance(v, float) and pd.isna(v):
                    r[k] = None
        return rows

    top_real = _top(ranked[~ranked["simulated"]])
    top_sim = _top(ranked[ranked["simulated"]])
    # Ranked order is kept within each phenotype, so head() is that phenotype's top.
    by_phen_real = {phen: _top(grp, TOP_N_PER_PHENOTYPE)
                    for phen, grp in ranked[~ranked["simulated"]].groupby("phenotype", sort=True)}
    by_phen_sim = {phen: _top(grp, TOP_N_PER_PHENOTYPE)
                   for phen, grp in ranked[ranked["simulated"]].groupby("phenotype", sort=True)}
    # Candidates whose evidence is a baseline difference between resilient and
    # susceptible groups, and significant: the design for a predictive marker.
    # Ranked like everything else (tier, then score); real studies only.
    real_ranked = ranked[~ranked["simulated"]]
    # Compared as text: a candidates.tsv round-trip gives "True"/"False", and a
    # missing value must read as False (bool(nan) is True).
    is_baseline = (real_ranked["is_predictive_baseline"].astype(str).str.strip().str.lower().eq("true")
                   if "is_predictive_baseline" in real_ranked.columns
                   else pd.Series(False, index=real_ranked.index))
    baseline = real_ranked[is_baseline]
    significant = baseline[pd.to_numeric(baseline["adjusted_p_value"], errors="coerce") <= SIGNIFICANCE_THRESHOLD]
    top_baseline = _top(significant, TOP_N_PREDICTIVE_BASELINE)
    predictive = {
        "n": int(len(baseline)),
        "n_significant": int(len(significant)),
        "n_significant_replicated": int((pd.to_numeric(significant["k_studies"], errors="coerce") >= 2).sum()),
        "studies": sorted({s for v in baseline["studies"].dropna() for s in str(v).split("|") if s}),
        "phenotypes": sorted(baseline["phenotype"].dropna().unique().tolist()),
        "top": top_baseline,
    }

    all_real = top_real + [r for rows in by_phen_real.values() for r in rows] + top_baseline
    all_rows = all_real + top_sim + [r for rows in by_phen_sim.values() for r in rows]
    for r in all_rows:
        r["card_page"] = _publish_card(r.pop("card_file", None), cards_dir)
    ann = _gene_annotations([r["feature_id_standardized"] for r in all_real
                             if r.get("feature_type") == "gene"])
    for r in all_real:
        r.update(ann.get(r["feature_id_standardized"], {}))

    return {
        "by_tier": by_tier,
        "by_evidence_class": by_class,
        "top_real": top_real,
        "top_simulated": top_sim,
        "top_real_by_phenotype": by_phen_real,
        "top_simulated_by_phenotype": by_phen_sim,
        "n_total": int(len(df)),
        "n_real": int((~df["simulated"]).sum()),
        "n_simulated": int(df["simulated"].sum()),
        "predictive_baseline": predictive,
    }


def _empty_predictive_baseline() -> dict:
    return {"n": 0, "n_significant": 0, "n_significant_replicated": 0, "studies": [], "phenotypes": [], "top": []}


def _taxid(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _species_label(taxid) -> str:
    term = species_term_for_taxid(taxid) if taxid is not None else None
    return term["scientific_name"] if term else ("unknown species" if taxid is None else f"taxid {taxid}")


def cross_species_gate(pools: list[dict]) -> list[dict]:
    """The WP6 dashboard gate, one row per phenotype and origin studied in two or more species.

    docs/mussel_clam_integration_plan.md, WP6: no cross-species pooling is shown
    until at least two species each have a within-species group with k >= 2 on
    the same phenotype. Before that, the orthogroup view is a lookup ("where is
    this oyster candidate's ortholog in mussels, and was it tested?"), not a
    meta-analysis. Only gene-level pools count; a species' pool reaches k >= 2
    when any of its genes was measured in two studies.
    """
    by_phenotype: dict[tuple, list[dict]] = {}
    for p in pools:
        if p["feature_type"] == "gene" and p.get("species_taxid") is not None:
            by_phenotype.setdefault((p["phenotype"], p["simulated"]), []).append(p)
    gate = []
    for (phen, sim), ps in sorted(by_phenotype.items()):
        if len({p["species_taxid"] for p in ps}) < 2:
            continue
        species = sorted(
            ({"species_taxid": p["species_taxid"], "species": p["species"], "max_k": p["max_k"],
              "studies": p["studies"]} for p in ps),
            key=lambda d: (-d["max_k"], d["species"]))
        n_k2 = sum(1 for d in species if d["max_k"] >= 2)
        is_open = n_k2 >= GATE_MIN_SPECIES_K2
        gate.append({
            "phenotype": phen, "simulated": bool(sim), "species": species,
            "n_species": len(species), "n_species_k2": n_k2, "open": is_open,
            "reason": (f"{n_k2} species each have a within-species pool with k >= 2"
                       if is_open else
                       f"only {n_k2} of {len(species)} species has a within-species pool with k >= 2; "
                       f"{GATE_MIN_SPECIES_K2} are needed"),
        })
    return gate


CROSS_SPECIES_COLUMNS = [
    "orthogroup_id", "phenotype", "feature_type", "copy_class", "tier", "evidence_class", "score",
    "species_taxids", "genes", "k_studies_by_species", "pooled_effect", "ci_lower", "ci_upper",
    "adjusted_p_value", "i_squared_between_species", "direction_consistency_species",
    "meets_replication_gate", "gates_failed", "simulated", "card_file",
]


def _cross_species_summary(gate: list[dict], phenotype_labels: dict, cards_dir: Path | None = None,
                           include_simulated: bool = False) -> dict:
    """Cross-species candidates, shown only for phenotypes past the gate.

    Rows for a phenotype whose gate is closed are counted as withheld and never
    ranked or linked: before the gate opens, a pooled cross-species estimate
    rests on single-study species and the page offers the ortholog lookup
    instead.
    """
    empty = {"available": CROSS_SPECIES_TABLE.exists(), "by_tier": [], "top_real": [],
             "top_simulated": [], "withheld": [], "n_shown": 0, "n_withheld": 0}
    if not CROSS_SPECIES_TABLE.exists():
        return empty
    df = pd.read_csv(CROSS_SPECIES_TABLE, sep="\t", low_memory=False, dtype={"species_taxids": str})
    df = _drop_simulated(df, include_simulated)
    if not len(df):
        return empty
    open_keys = {(g["phenotype"], g["simulated"]) for g in gate if g["open"]}
    keys = list(zip(df["phenotype"], df["simulated"]))
    shown_mask = pd.Series([k in open_keys for k in keys], index=df.index)
    shown, held = df[shown_mask], df[~shown_mask]

    from prioritize.cross_species import TIER_ORDER as CROSS_TIERS

    rank = {t: i for i, t in enumerate(CROSS_TIERS)}
    shown = (shown.assign(_r=shown["tier"].map(rank).fillna(len(CROSS_TIERS)))
             .sort_values(["_r", "score"], ascending=[True, False]))
    cols = [c for c in CROSS_SPECIES_COLUMNS if c in shown.columns]

    def _top(frame: pd.DataFrame) -> list[dict]:
        rows = frame.head(TOP_N_CROSS_SPECIES)[cols].to_dict(orient="records")
        for r in rows:
            r["phenotype_label"] = phenotype_labels.get(r["phenotype"], r["phenotype"])
            k_by = dict(pair.split(":", 1) for pair in str(r.get("k_studies_by_species") or "").split("|")
                        if ":" in pair)
            r["species_genes"] = [
                {"species_taxid": _taxid(t), "species": _species_label(_taxid(t)), "gene": g,
                 "k_studies": _taxid(k_by.get(t))}
                for t, g in (pair.split(":", 1) for pair in str(r.get("genes") or "").split("|") if ":" in pair)]
            for k, v in list(r.items()):
                if isinstance(v, float) and pd.isna(v):
                    r[k] = None
            r["card_page"] = _publish_card(r.pop("card_file", None), cards_dir)
        return rows

    return {
        "available": True,
        "by_tier": [{"simulated": bool(s), "phenotype": ph, "tier": t, "n": int(n)}
                    for (s, ph, t), n in shown.groupby(["simulated", "phenotype", "tier"]).size().items()],
        "top_real": _top(shown[~shown["simulated"]]),
        "top_simulated": _top(shown[shown["simulated"]]),
        "withheld": [{"simulated": bool(s), "phenotype": ph, "n": int(n)}
                     for (s, ph), n in held.groupby(["simulated", "phenotype"]).size().items()],
        "n_shown": int(len(shown)),
        "n_withheld": int(len(held)),
    }


def _attach_ortholog_lookup(rows: list[dict], evidence: pd.DataFrame) -> None:
    """For each single-species candidate, the other species' genes in its orthogroup.

    This is the lookup the gate always allows: where a candidate's ortholog sits
    in another species, and whether that species tested it, under which
    phenotypes. No estimate is pooled here.
    """
    if not rows or not len(evidence) or "orthogroup_id" not in evidence.columns:
        for r in rows:
            r["orthologs"] = []
        return
    wanted = {str(r.get("orthogroup_id")) for r in rows if r.get("orthogroup_id")}
    ev = evidence[evidence["orthogroup_id"].astype(str).isin(wanted)]
    adj = pd.to_numeric(ev.get("adjusted_p_value"), errors="coerce")
    ev = ev.assign(_taxid=ev["species_taxid"].map(_taxid), _sig=adj <= SIGNIFICANCE_THRESHOLD)
    index: dict[tuple, list[dict]] = {}
    for (og, sim, taxid), grp in ev.groupby(["orthogroup_id", "simulated", "_taxid"], sort=True):
        index.setdefault((str(og), bool(sim)), []).append({
            "species_taxid": taxid, "species": _species_label(taxid),
            "genes": sorted({str(g) for g in grp["feature_id_standardized"].dropna()}),
            "phenotypes": sorted({str(p) for p in grp["phenotype"].dropna()}),
            "studies": sorted({str(s) for s in grp["study_id"].dropna()}),
            "n_significant_records": int(grp["_sig"].sum()),
        })
    for r in rows:
        own = r.get("species_taxid")
        r["orthologs"] = [o for o in index.get((str(r.get("orthogroup_id")), bool(r.get("simulated"))), [])
                          if o["species_taxid"] != own]


def _publish_card(card_file, cards_dir: Path | None) -> str | None:
    """Copy one evidence card next to the dashboard so the site can link to it.

    Returns the page path relative to the home page (index.qmd), which loads
    data.json from one directory below it, so the link includes that directory
    (e.g. ``dashboard/cards/<card>.html``). None when the candidate has no card.
    """
    if not card_file or cards_dir is None or (isinstance(card_file, float) and pd.isna(card_file)):
        return None
    src = Path(str(card_file))
    if not src.is_absolute():
        src = REPO_ROOT / src
    if not src.exists():
        return None
    cards_dir.mkdir(parents=True, exist_ok=True)
    dest = cards_dir / src.name
    lines = src.read_text().splitlines()
    if lines and lines[0].startswith("# "):
        # Promote the card's H1 to Quarto front matter so the site does not show it twice.
        title = lines[0][2:].strip().replace('"', "'")
        body = "\n".join(lines[1:]).lstrip("\n")
    else:
        title = src.stem
        body = "\n".join(lines)
    # Cards are reached from the home page, which has no sidebar; match it.
    dest.write_text(
        f'---\ntitle: "{title}"\nsidebar: false\nbread-crumbs: false\n---\n\n{body}\n'
    )
    return f"{cards_dir.parent.name}/{cards_dir.name}/{src.stem}.html"


def _pool_summary(include_simulated: bool = False) -> list[dict]:
    """One row per (phenotype, feature_type, origin, species) pool actually formed.

    Gene-level pools never cross species, so a phenotype studied in two species
    is two pools, each with its own study count. The WP6 gate reads them.
    """
    pools = []
    if not META_ANALYSIS_DIR.exists():
        return pools
    seen = set()
    # A per-phenotype run and the all-phenotypes run describe the same pools;
    # read the specific files first so `seen` credits each pool to one of them.
    # Orthogroup-level (cross-species) results share the suffix but not the
    # columns; they are summarized separately, not as gene-level pools.
    paths = sorted((p for p in META_ANALYSIS_DIR.glob("*_meta_analysis.tsv")
                    if not p.name.endswith("_orthogroup_meta_analysis.tsv")),
                   key=lambda p: (p.name.startswith("all_phenotypes"), p.name))
    for path in paths:
        df = pd.read_csv(path, sep="\t", low_memory=False)
        if not len(df):
            continue
        df = _drop_simulated(df, include_simulated)
        if not len(df):
            continue
        adj = pd.to_numeric(df["adjusted_p_value"], errors="coerce")
        k = pd.to_numeric(df["k_studies"], errors="coerce").fillna(0)
        if "species_taxid" not in df.columns:
            df["species_taxid"] = None
        df["_taxid"] = df["species_taxid"].map(_taxid)
        # A baseline contrast and an exposure response under one phenotype are
        # separate pools (meta_analysis.run.GROUP_KEYS).
        if "contrast_type" not in df.columns:
            df["contrast_type"] = DEFAULT_CONTRAST_TYPE
        for (phen, ftype, sim, taxid, contrast), grp in df.groupby(
                ["phenotype", "feature_type", "simulated", "_taxid", "contrast_type"], dropna=False):
            key = (phen, ftype, bool(sim), taxid, contrast)
            if key in seen:
                continue
            seen.add(key)
            idx = grp.index
            studies: set[str] = set()
            for s in grp["studies"].dropna():
                studies.update(str(s).split("|"))
            pools.append({
                "phenotype": phen,
                "feature_type": ftype,
                "simulated": bool(sim),
                "species_taxid": taxid,
                "contrast_type": contrast,
                "species": _species_label(taxid),
                "studies": sorted(studies),
                "n_features": int(len(grp)),
                "n_features_k2plus": int((k.loc[idx] >= 2).sum()),
                "n_features_k2plus_significant": int(
                    ((k.loc[idx] >= 2) & (adj.loc[idx] <= SIGNIFICANCE_THRESHOLD)).sum()
                ),
                "max_k": int(k.loc[idx].max()),
                "source_file": _repo_relative(path),
            })
    return pools


# --------------------------------------------------------------------------- #
# Repository activity, tests, implementation status
# --------------------------------------------------------------------------- #


def _activity() -> dict:
    log = _git("log", "--format=%H%x09%ad%x09%s", "--date=short")
    commits = []
    for line in log.splitlines():
        parts = line.split("\t", 2)
        if len(parts) == 3:
            commits.append({"sha": parts[0][:7], "date": parts[1], "subject": parts[2]})
    by_week: Counter = Counter()
    for c in commits:
        d = _dt.date.fromisoformat(c["date"])
        week_start = d - _dt.timedelta(days=d.weekday())
        by_week[week_start.isoformat()] += 1
    # Fill empty weeks so the activity chart has a continuous axis.
    weeks = []
    if by_week:
        start = _dt.date.fromisoformat(min(by_week))
        end = _dt.date.fromisoformat(max(by_week))
        cur = start
        while cur <= end:
            weeks.append({"week": cur.isoformat(), "n_commits": by_week.get(cur.isoformat(), 0)})
            cur += _dt.timedelta(days=7)
    return {
        "n_commits": len(commits),
        "first_commit_date": commits[-1]["date"] if commits else None,
        "last_commit_date": commits[0]["date"] if commits else None,
        "commits_by_week": weeks,
        "recent_commits": commits[:12],
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD") or None,
    }


def _test_summary() -> dict:
    files = sorted(TESTS_DIR.glob("test_*.py"))
    n_tests = 0
    for f in files:
        n_tests += len(re.findall(r"^\s*def test_", f.read_text(), flags=re.M))
    return {"n_test_files": len(files), "n_tests": n_tests}


def _implementation_status() -> dict:
    """Count table rows under each section of the candid status page."""
    if not IMPLEMENTATION_STATUS_MD.exists():
        return {"sections": [], "source": None}
    sections = []
    current = None
    for line in IMPLEMENTATION_STATUS_MD.read_text().splitlines():
        if line.startswith("## "):
            current = {"title": line[3:].strip(), "n_items": 0}
            sections.append(current)
        elif current is not None:
            stripped = line.strip()
            is_table_row = (stripped.startswith("|") and not stripped.startswith("|---")
                            and not re.match(r"^\|\s*Component\s*\|", stripped))
            is_bullet = stripped.startswith("- ")
            if is_table_row or is_bullet:
                current["n_items"] += 1
    return {"sections": sections, "source": "docs/implementation_status.md"}


# --------------------------------------------------------------------------- #
# Assembly
# --------------------------------------------------------------------------- #


def build_dashboard_data(out_path: Path | None = None, include_simulated: bool = False) -> dict:
    """Summarize the registry, pipeline and findings into the dashboard's data.json.

    Simulated demo studies and their evidence are left out unless
    include_simulated is set: they are pipeline fixtures, not findings.
    """
    out_path = Path(out_path) if out_path else DEFAULT_OUT_PATH
    phenotype_vocab = load_vocab("phenotype_ontology")
    phenotype_labels = {k: v["label"] for k, v in phenotype_vocab.items()}
    phenotype_relevance = {k: v.get("resilience_relevance") for k, v in phenotype_vocab.items()}
    stressor_labels = {k: v["label"] for k, v in load_vocab("stressor_ontology").items()}

    studies = _load_studies(include_simulated)
    evidence = _load_evidence(include_simulated)
    study_rows, comparison_rows = _study_rows(studies, evidence)

    real_studies = [s for s in study_rows if not s["simulated"]]
    sim_studies = [s for s in study_rows if s["simulated"]]
    real_comparisons = [c for c in comparison_rows if not c["simulated"]]

    cards_dir = out_path.parent / CARDS_SUBDIR
    if cards_dir.exists():
        for stale in cards_dir.glob("*.md"):
            stale.unlink()
    candidates = _candidate_summary(phenotype_labels, cards_dir=cards_dir, include_simulated=include_simulated)
    pools = _pool_summary(include_simulated)
    gate = cross_species_gate(pools)
    cross_species = _cross_species_summary(gate, phenotype_labels, cards_dir=cards_dir,
                                           include_simulated=include_simulated)
    cross_species["gate"] = gate
    seen_rows: dict[int, dict] = {}
    for key in ("top_real", "top_simulated"):
        for r in candidates[key]:
            seen_rows[id(r)] = r
    for rows in [*candidates["top_real_by_phenotype"].values(), *candidates["top_simulated_by_phenotype"].values()]:
        for r in rows:
            seen_rows[id(r)] = r
    _attach_ortholog_lookup(list(seen_rows.values()), evidence)

    def _tier_count(simulated: bool, tier: str) -> int:
        return sum(r["n"] for r in candidates["by_tier"]
                   if r["simulated"] == simulated and r["tier"] == tier)

    real_resilience_comparisons = [
        c for c in real_comparisons if c["resilience_classification"] == "resilience"
    ]
    real_stressors = sorted({c["stressor"] for c in real_comparisons if c["stressor"]})

    kpis = {
        "n_studies_real": len(real_studies),
        "n_studies_simulated": len(sim_studies),
        "n_comparisons_real": len(real_comparisons),
        "n_comparisons_real_harmonized": sum(c["harmonized"] for c in real_comparisons),
        "n_comparisons_simulated": len(comparison_rows) - len(real_comparisons),
        "n_evidence_records_real": int((~evidence["simulated"]).sum()) if len(evidence) else 0,
        "n_evidence_records_simulated": int(evidence["simulated"].sum()) if len(evidence) else 0,
        "n_candidates_real_high_priority": _tier_count(False, "high_priority_cross_study"),
        "n_candidates_real_multi_omics": _tier_count(False, "multi_omics_convergence"),
        "n_candidates_real_metabolite_linked": _tier_count(False, "metabolite_linked_convergence"),
        "n_candidates_real_emerging": _tier_count(False, "emerging"),
        "n_real_stressor_classes": len(real_stressors),
        "real_stressor_classes": real_stressors,
        "n_real_comparisons_with_resilience_phenotype": len(real_resilience_comparisons),
        "n_real_pools_k2plus": sum(1 for p in pools if not p["simulated"] and p["max_k"] >= 2),
        "n_real_species_with_evidence": len({p["species_taxid"] for p in pools if not p["simulated"]}),
        "n_real_phenotypes_cross_species_open": sum(1 for g in gate if not g["simulated"] and g["open"]),
    }

    head_sha = _git("rev-parse", "HEAD")
    data = {
        "generated_at": _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat(),
        "repo_url": REPO_URL,
        "commit": {
            "sha": head_sha or None,
            "short": head_sha[:7] if head_sha else None,
            "date": _git("log", "-1", "--format=%ad", "--date=short") or None,
            "subject": _git("log", "-1", "--format=%s") or None,
        },
        "artifacts_present": {
            "evidence_table": EVIDENCE_TABLE.exists(),
            "meta_analysis": META_ANALYSIS_DIR.exists() and any(META_ANALYSIS_DIR.glob("*.tsv")),
            "candidates": CANDIDATES_TABLE.exists(),
            "cross_species_candidates": CROSS_SPECIES_TABLE.exists(),
            "top_candidates_summary": TOP_CANDIDATES_MD.exists(),
            "real_crosswalk": bool(_built_crosswalks()),
        },
        "crosswalks": [
            {
                "species": entry["species"],
                "ncbi_taxid": entry["ncbi_taxid"],
                "path": str(entry["path"].relative_to(REPO_ROOT))
                if entry["path"].is_relative_to(REPO_ROOT) else str(entry["path"]),
                "assemblies": entry["assemblies"],
                "built": entry["path"].exists(),
            }
            for entry in registered_crosswalks()
        ],
        "includes_simulated": include_simulated,
        "kpis": kpis,
        "studies": study_rows,
        "comparisons": comparison_rows,
        "evidence": _evidence_summary(evidence),
        "candidates": candidates,
        "pools": pools,
        "cross_species": cross_species,
        "species": species_summary(studies, evidence),
        "activity": _activity(),
        "tests": _test_summary(),
        "implementation_status": _implementation_status(),
        "vocab": {
            "phenotype_labels": phenotype_labels,
            "phenotype_resilience_relevance": phenotype_relevance,
            "stressor_labels": stressor_labels,
        },
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(data, fh, indent=1, default=_json_default)
    return data


def _json_default(obj):
    if isinstance(obj, (pd.Timestamp, _dt.date, _dt.datetime)):
        return obj.isoformat()
    if hasattr(obj, "item"):
        return obj.item()
    if isinstance(obj, set):
        return sorted(obj)
    return str(obj)
