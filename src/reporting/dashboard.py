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

from common import REPO_ROOT, REPORTS_DIR, STUDIES_DIR, load_vocab

EVIDENCE_TABLE = REPORTS_DIR / "evidence" / "evidence_table.tsv"
CANDIDATES_TABLE = REPORTS_DIR / "evidence_cards" / "candidates.tsv"
META_ANALYSIS_DIR = REPORTS_DIR / "meta_analysis"
MANIFESTS_DIR = REPORTS_DIR / "manifests"
TOP_CANDIDATES_MD = REPORTS_DIR / "top_candidates_summary.md"
IMPLEMENTATION_STATUS_MD = REPO_ROOT / "docs" / "implementation_status.md"
REAL_CROSSWALK = REPO_ROOT / "data" / "reference" / "crosswalk" / "mgigas_gene_id_crosswalk.tsv"
TESTS_DIR = REPO_ROOT / "tests"
DEFAULT_OUT_PATH = REPO_ROOT / "docs" / "dashboard" / "data.json"
CARDS_SUBDIR = "cards"  # evidence cards for the top candidates, copied beside data.json
REPO_URL = "https://github.com/sr320/AREE"

SIGNIFICANCE_THRESHOLD = 0.05
TOP_N_CANDIDATES = 15


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


def _load_studies() -> list[dict]:
    studies = []
    for path in sorted(STUDIES_DIR.glob("*.yaml")):
        if path.name.startswith("_"):
            continue
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


def _load_evidence() -> pd.DataFrame:
    if not EVIDENCE_TABLE.exists():
        return pd.DataFrame()
    df = pd.read_csv(EVIDENCE_TABLE, sep="\t", low_memory=False)
    df["simulated"] = df["simulated"].map(_bool)
    return df


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


def _gene_annotations(gene_ids: list[str]) -> dict[str, dict]:
    """Symbol and description for real NCBI GeneIDs, from the committed crosswalk."""
    if not REAL_CROSSWALK.exists() or not gene_ids:
        return {}
    wanted = set(map(str, gene_ids))
    out: dict[str, dict] = {}
    for chunk in pd.read_csv(REAL_CROSSWALK, sep="\t", comment="#", dtype=str,
                             chunksize=20000):
        hit = chunk[chunk["ncbi_gene_id"].isin(wanted)]
        for _, row in hit.iterrows():
            out[row["ncbi_gene_id"]] = {
                "gene_symbol": row.get("gene_symbol") or "",
                "gene_description": row.get("gene_description") or "",
                "gene_type": row.get("gene_type") or "",
            }
        if len(out) >= len(wanted):
            break
    return out


def _candidate_summary(phenotype_labels: dict, cards_dir: Path | None = None) -> dict:
    if not CANDIDATES_TABLE.exists():
        return {"by_tier": [], "by_evidence_class": [], "top_real": [], "top_simulated": [],
                "n_total": 0}
    df = pd.read_csv(CANDIDATES_TABLE, sep="\t", low_memory=False)
    df["simulated"] = df["simulated"].map(_bool)
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
            "context_replication", "score", "k_studies", "studies", "total_sample_size",
            "pooled_effect", "adjusted_p_value", "direction_consistency", "i_squared",
            "n_supporting_layers", "mapping_confidences", "simulated", "card_file"]
    cols = [c for c in cols if c in ranked.columns]

    def _top(frame: pd.DataFrame) -> list[dict]:
        rows = frame.head(TOP_N_CANDIDATES)[cols].to_dict(orient="records")
        for r in rows:
            r["feature_id_standardized"] = str(r["feature_id_standardized"])
            r["phenotype_label"] = phenotype_labels.get(r.get("phenotype"), r.get("phenotype"))
            for k, v in list(r.items()):
                if isinstance(v, float) and pd.isna(v):
                    r[k] = None
        return rows

    top_real = _top(ranked[~ranked["simulated"]])
    top_sim = _top(ranked[ranked["simulated"]])
    for r in top_real + top_sim:
        r["card_page"] = _publish_card(r.pop("card_file", None), cards_dir)
    ann = _gene_annotations([r["feature_id_standardized"] for r in top_real
                             if r.get("feature_type") == "gene"])
    for r in top_real:
        r.update(ann.get(r["feature_id_standardized"], {}))

    return {
        "by_tier": by_tier,
        "by_evidence_class": by_class,
        "top_real": top_real,
        "top_simulated": top_sim,
        "n_total": int(len(df)),
        "n_real": int((~df["simulated"]).sum()),
        "n_simulated": int(df["simulated"].sum()),
    }


def _publish_card(card_file, cards_dir: Path | None) -> str | None:
    """Copy one evidence card next to the dashboard so the site can link to it.

    Returns the site-relative page path (without extension, as Quarto links it),
    or None when the candidate has no rendered card.
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
    return f"{cards_dir.name}/{src.stem}.html"


def _pool_summary() -> list[dict]:
    """One row per (phenotype, feature_type, origin) pool actually formed."""
    pools = []
    if not META_ANALYSIS_DIR.exists():
        return pools
    seen = set()
    # A per-phenotype run and the all-phenotypes run describe the same pools;
    # read the specific files first so `seen` credits each pool to one of them.
    paths = sorted(META_ANALYSIS_DIR.glob("*_meta_analysis.tsv"),
                   key=lambda p: (p.name.startswith("all_phenotypes"), p.name))
    for path in paths:
        df = pd.read_csv(path, sep="\t", low_memory=False)
        if not len(df):
            continue
        df["simulated"] = df["simulated"].map(_bool)
        adj = pd.to_numeric(df["adjusted_p_value"], errors="coerce")
        k = pd.to_numeric(df["k_studies"], errors="coerce").fillna(0)
        for (phen, ftype, sim), grp in df.groupby(["phenotype", "feature_type", "simulated"]):
            key = (phen, ftype, bool(sim))
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


def build_dashboard_data(out_path: Path | None = None) -> dict:
    out_path = Path(out_path) if out_path else DEFAULT_OUT_PATH
    phenotype_vocab = load_vocab("phenotype_ontology")
    phenotype_labels = {k: v["label"] for k, v in phenotype_vocab.items()}
    phenotype_relevance = {k: v.get("resilience_relevance") for k, v in phenotype_vocab.items()}
    stressor_labels = {k: v["label"] for k, v in load_vocab("stressor_ontology").items()}

    studies = _load_studies()
    evidence = _load_evidence()
    study_rows, comparison_rows = _study_rows(studies, evidence)

    real_studies = [s for s in study_rows if not s["simulated"]]
    sim_studies = [s for s in study_rows if s["simulated"]]
    real_comparisons = [c for c in comparison_rows if not c["simulated"]]

    cards_dir = out_path.parent / CARDS_SUBDIR
    if cards_dir.exists():
        for stale in cards_dir.glob("*.md"):
            stale.unlink()
    candidates = _candidate_summary(phenotype_labels, cards_dir=cards_dir)
    pools = _pool_summary()

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
        "n_candidates_real_emerging": _tier_count(False, "emerging"),
        "n_real_stressor_classes": len(real_stressors),
        "real_stressor_classes": real_stressors,
        "n_real_comparisons_with_resilience_phenotype": len(real_resilience_comparisons),
        "n_real_pools_k2plus": sum(1 for p in pools if not p["simulated"] and p["max_k"] >= 2),
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
            "top_candidates_summary": TOP_CANDIDATES_MD.exists(),
            "real_crosswalk": REAL_CROSSWALK.exists(),
        },
        "kpis": kpis,
        "studies": study_rows,
        "comparisons": comparison_rows,
        "evidence": _evidence_summary(evidence),
        "candidates": candidates,
        "pools": pools,
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
