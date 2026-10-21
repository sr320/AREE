"""The read-only query layer, the site export, and the MCP server over the demo.

Every test reads a temporary demo run (isolated_reports), never reports/.
"""
from __future__ import annotations

import asyncio
import json
from pathlib import Path

import pandas as pd
import pytest
from conftest import ALL_DEMO_STUDY_IDS, harmonize_all_demo_studies

import query
from common import study_files
from reporting.evidence_cards import build_evidence_cards

CONFLICT_GENE = "LOC105331241"  # sod1 in the demo: opposite directions across larval studies


@pytest.fixture
def sources(isolated_reports):
    return query.Sources(reports_dir=isolated_reports["evidence_table_path"].parents[1])


@pytest.fixture
def demo_sources(sources):
    harmonize_all_demo_studies()
    build_evidence_cards(phenotype="larval_viability", feature_type="gene")
    return sources


def test_missing_outputs_are_reported_not_hidden(sources):
    out = query.overview(sources)
    assert out["evidence_records"] == {"real": 0, "simulated": 0}
    assert out["artifacts"]["evidence_table"]["present"] is False
    assert "aree harmonize" in out["artifacts"]["evidence_table"]["produced_by"]
    assert query.get_evidence(feature="x", sources=sources)["artifacts"]["evidence_table"]["present"] is False
    assert query.list_candidates(sources=sources)["n_candidates"] == 0


def test_studies_come_from_every_yaml_and_filter(sources):
    every = query.list_studies(sources=sources)
    # Real (registry/studies/) and demo (registry/demo_studies/) alike.
    assert every["n_studies"] == len(study_files())
    real_only = query.list_studies(include_simulated=False, sources=sources)
    assert not {s["study_id"] for s in real_only["studies"]} & set(ALL_DEMO_STUDY_IDS)
    heat = query.list_studies(stressor="temperature", sources=sources)
    assert heat["n_studies"] > 0
    assert all("temperature" in s["stressors"] for s in heat["studies"])


def test_get_study_returns_the_record_or_known_ids(sources):
    study_id = ALL_DEMO_STUDY_IDS[0]
    out = query.get_study(study_id, sources)
    assert out["study"]["study_id"] == study_id
    assert out["registration_file"].endswith(f"{study_id}.yaml")
    missing = query.get_study("NOT_A_STUDY", sources)
    assert "error" in missing and study_id in missing["known_study_ids"]


def test_evidence_keeps_provenance_and_surfaces_conflicts(demo_sources):
    out = query.get_evidence(feature=CONFLICT_GENE, phenotype="larval_viability", sources=demo_sources)
    assert out["n_records"] >= 2
    assert out["conflicting_directions"] is True
    for rec in out["records"]:
        for field in ("evidence_id", "study_id", "simulated", "source_file", "workflow_version",
                      "date_generated", "mapping_confidence"):
            assert field in rec
        assert rec["simulated"] is True or str(rec["simulated"]).lower() == "true"
    json.dumps(out)  # every value serializes: no NaN, no numpy scalars
    assert "not a validated biomarker" in out["note"]


def test_replicated_features_are_counted_by_study_with_their_directions(demo_sources):
    out = query.find_replicated_features(phenotype="larval_viability", max_adjusted_p=1.0, limit=500,
                                         sources=demo_sources)
    hit = next(f for f in out["features"] if f["feature_id_standardized"] == CONFLICT_GENE)
    assert hit["n_studies"] >= 2 and len(hit["studies"]) == hit["n_studies"]
    assert hit["n_studies_up"] >= 1 and hit["n_studies_down"] >= 1 and hit["conflicting_directions"]
    assert out["n_conflicting"] >= 1 and all(f["n_studies"] >= 2 for f in out["features"])
    # Most studies first, and among equals, consistent before conflicting.
    keys = [(-f["n_studies"], f["conflicting_directions"]) for f in out["features"]]
    assert keys == sorted(keys)
    assert "not fully independent" in out["note"]
    json.dumps(out)

    stricter = query.find_replicated_features(phenotype="larval_viability", min_studies=50,
                                              sources=demo_sources)
    assert stricter["n_features"] == 0 and stricter["features"] == []
    real_only = query.find_replicated_features(max_adjusted_p=1.0, include_simulated=False,
                                               sources=demo_sources)
    assert real_only["n_features"] == 0


def test_simulated_evidence_can_be_excluded(demo_sources):
    assert query.get_evidence(feature=CONFLICT_GENE, include_simulated=False, sources=demo_sources)["n_records"] == 0


def test_search_features_summarizes_one_row_per_feature(demo_sources):
    out = query.search_features(CONFLICT_GENE.lower(), sources=demo_sources)
    hit = next(f for f in out["features"] if f["feature_id_standardized"] == CONFLICT_GENE)
    assert hit["n_studies"] >= 2 and hit["any_simulated"]
    assert {"up", "down"} <= set(hit["directions"])
    json.dumps(out)


def test_candidates_and_their_evidence_card(demo_sources):
    cands = query.list_candidates(phenotype="larval_viability", sources=demo_sources)
    assert cands["n_candidates"] > 0 and sum(cands["tier_counts"].values()) == cands["n_candidates"]
    out = query.get_evidence_card(CONFLICT_GENE, phenotype="larval_viability", sources=demo_sources)
    assert out["n_cards"] >= 1
    assert any(c["markdown"] and "Conflicting direction" in c["markdown"] for c in out["cards"])
    json.dumps(out)
    assert query.get_evidence_card("NO_SUCH_GENE", sources=demo_sources)["n_cards"] == 0


def test_vocabulary_lookup():
    out = query.list_vocabulary("phenotype")
    assert {"survival", "thermal_tolerance"} <= {t["id"] for t in out["terms"]}
    assert "error" in query.list_vocabulary("nonsense")


def test_site_export_writes_sorted_parquet_and_meta(demo_sources, tmp_path):
    pytest.importorskip("pyarrow")
    from query.site_export import export_site_data

    out_dir = tmp_path / "site"
    meta = export_site_data(out_dir, sources=demo_sources, include_simulated=True)
    evidence = pd.read_parquet(out_dir / "evidence.parquet")
    assert len(evidence) == meta["tables"]["evidence"]["total"] > 0
    assert evidence["simulated"].dtype == bool
    ids = evidence["feature_id_standardized"].dropna().tolist()
    assert ids == sorted(ids)  # sorted, so a gene lookup touches few row groups
    assert set(query.EVIDENCE_COLUMNS) <= set(evidence.columns)
    assert all(isinstance(v, str) for v in evidence["feature_id_standardized"].dropna())  # ids stay text
    # The feature index resolves every evidence row's identifiers, with counts that add up.
    features = pd.read_parquet(out_dir / "features.parquet")
    assert features["n_records"].sum() == len(evidence)
    assert CONFLICT_GENE in set(features["feature_id_standardized"])
    assert not (features[["feature_id_standardized", "feature_id_original"]] == "nan").any().any()
    candidates = pd.read_parquet(out_dir / "candidates.parquet")
    assert len(candidates) == meta["tables"]["candidates"]["total"]
    # Real booleans, so `WHERE is_predictive_baseline` works in the page's SQL.
    assert candidates["is_predictive_baseline"].dtype == bool
    assert {"contrast_type", "gene_symbol", "gene_description"} <= set(candidates.columns)
    assert meta["facets"]["candidates"]["contrast_type"] == ["exposure_response"]
    assert len(pd.read_parquet(out_dir / "studies.parquet")) == meta["tables"]["studies"]["total"]
    on_disk = json.loads((out_dir / "meta.json").read_text())
    assert "larval_viability" in on_disk["facets"]["evidence"]["phenotype"]
    assert on_disk["note"] == query.INTERPRETATION_NOTE

    # The published site leaves simulated demo data out unless asked.
    real_only = export_site_data(tmp_path / "real_only", sources=demo_sources)
    assert real_only["tables"]["evidence"]["total"] == 0
    assert real_only["tables"]["studies"]["simulated"] == 0
    assert real_only["tables"]["studies"]["real"] > 0


def test_mcp_server_exposes_the_tools_and_answers(demo_sources):
    pytest.importorskip("mcp")
    from query.mcp_server import build_server

    server = build_server(demo_sources)

    async def exercise():
        tools = {t.name for t in await server.list_tools()}
        result = await server.call_tool("get_evidence", {"feature": CONFLICT_GENE, "phenotype": "larval_viability"})
        return tools, result

    tools, result = asyncio.run(exercise())
    assert tools == {"aree_overview", "list_vocabulary", "list_studies", "get_study", "search_features",
                     "get_evidence", "find_replicated_features", "list_candidates", "get_evidence_card"}
    # Depending on the mcp 1.x release, call_tool returns the content blocks or (content, structured).
    content = result[0] if isinstance(result, tuple) else result
    payload = json.loads(content[0].text)
    assert payload["conflicting_directions"] is True


def _comparable(out):
    """An answer with the parts that legitimately differ between sources removed:
    where outputs live, and card paths (repo paths vs file names)."""
    if isinstance(out, dict):
        return {k: (Path(v).name if k == "card_file" and v else _comparable(v))
                for k, v in out.items() if k != "artifacts"}
    if isinstance(out, list):
        return [_comparable(v) for v in out]
    return out


def _by_evidence_id(out):
    out["records"] = sorted(out["records"], key=lambda r: r["evidence_id"])
    return out


@pytest.fixture
def site_sources(demo_sources, tmp_path):
    pytest.importorskip("pyarrow")
    from query.site_export import export_site_data

    export_site_data(tmp_path / "site", sources=demo_sources, include_simulated=True)
    return query.Sources.from_site_data(tmp_path / "site")


def test_published_tables_answer_exactly_as_the_reports_do(demo_sources, site_sources):
    # Evidence is re-sorted for export, so records tied on adjusted p are compared as sets.
    questions = [
        lambda s: query.overview(s),
        lambda s: query.search_features("LOC", limit=500, sources=s),
        lambda s: query.search_features(CONFLICT_GENE.lower(), include_simulated=False, sources=s),
        lambda s: _by_evidence_id(query.get_evidence(feature=CONFLICT_GENE, sources=s)),
        lambda s: _by_evidence_id(query.get_evidence(phenotype="larval_viability", max_adjusted_p=0.05,
                                                     species="gigas", limit=500, sources=s)),
        lambda s: query.list_candidates(limit=500, sources=s),
        lambda s: query.find_replicated_features(max_adjusted_p=1.0, limit=500, sources=s),
        lambda s: query.list_candidates(phenotype="larval_viability", tier="emerging", sources=s),
        lambda s: query.get_evidence_card(CONFLICT_GENE, phenotype="larval_viability", sources=s),
        lambda s: query.get_study(ALL_DEMO_STUDY_IDS[0], s)["n_evidence_records"],
    ]
    for ask in questions:
        from_reports, from_site = ask(demo_sources), ask(site_sources)
        assert _comparable(from_reports) == _comparable(from_site)
        json.dumps(from_site)
    card = query.get_evidence_card(CONFLICT_GENE, phenotype="larval_viability", sources=site_sources)
    assert card["cards"][0]["markdown"]  # served from cards.parquet, not reports/


def test_published_tables_are_queried_safely_from_many_threads(site_sources):
    from concurrent.futures import ThreadPoolExecutor

    def ask(_):
        return query.get_evidence(feature=CONFLICT_GENE, sources=site_sources)["n_records"]

    with ThreadPoolExecutor(max_workers=8) as pool:
        counts = set(pool.map(ask, range(32)))
    assert len(counts) == 1 and counts.pop() >= 2


def test_fetch_site_data_copies_every_published_table(site_sources, tmp_path):
    from query.site_export import fetch_site_data

    published = site_sources.data_dir
    meta = fetch_site_data(tmp_path / "fetched", base_url=published.as_uri())
    for name in [t["file"] for t in meta["tables"].values()] + ["cards.parquet", "meta.json"]:
        assert (tmp_path / "fetched" / name).read_bytes() == (published / name).read_bytes()
    fetched = query.Sources.from_site_data(tmp_path / "fetched")
    assert query.overview(fetched)["evidence_records"] == query.overview(site_sources)["evidence_records"]


def test_opposing_methylation_directions_count_as_conflicting():
    assert query._conflicting({"hyper": 3, "hypo": 1})
    assert query._conflicting({"up": 1, "down": 1})
    assert not query._conflicting({"up": 2, "hyper": 1})
