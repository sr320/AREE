"""The Ask AREE chat: its tools, its turn loop, its limits and its web API.

Claude is replaced by a scripted fake client, so these tests never call the API.
The data are a few synthetic records written as published Parquet tables, marked
real (simulated=False) because the chat serves real studies only, plus one
simulated record that must never reach a visitor.
"""
from __future__ import annotations

import copy
import json
from types import SimpleNamespace

import pandas as pd
import pytest

pytest.importorskip("pyarrow")
pytest.importorskip("fastapi")
pytest.importorskip("anthropic")

from conftest import ALL_DEMO_STUDY_IDS  # noqa: E402
from fastapi.testclient import TestClient  # noqa: E402

import query  # noqa: E402
from chat import tools  # noqa: E402
from chat.conversation import ChatConfig, Conversation, Conversations, run_turn  # noqa: E402
from chat.limits import DailyBudget, RateLimiter, cost_usd  # noqa: E402
from chat.server import Settings, create_app  # noqa: E402

GENE = "LOC900001"


def _evidence_row(evidence_id, study_id, feature, direction, adj_p, phenotype="thermal_tolerance",
                  simulated=False, feature_type="gene"):
    return {
        "evidence_id": evidence_id, "study_id": study_id, "comparison_id": "c1",
        "contrast_type": "exposure_response", "simulated": simulated,
        "species": "Crassostrea gigas", "feature_id_original": feature,
        "feature_id_standardized": feature, "feature_type": feature_type, "orthogroup_id": None,
        "mapping_confidence": "exact", "molecular_direction": direction, "effect_size": 1.2,
        "effect_size_type": "log2FC", "standard_error": 0.3, "ci_lower": 0.6, "ci_upper": 1.8,
        "p_value": adj_p / 10, "adjusted_p_value": adj_p, "sample_size": 12, "tissue": "gill",
        "life_stage": "adult", "stressor": "temperature", "phenotype": phenotype,
        "phenotype_direction": "higher_is_better", "resilience_classification": "resilience_associated",
        "analysis_method": "DESeq2", "quality_flags": None, "genome_assembly": "GCF_963853765.1",
        "annotation_version": "RS_2024", "source_file": "results.tsv", "workflow_version": "0.1.0",
        "date_generated": "2026-01-01",
    }


@pytest.fixture
def site(tmp_path):
    """Published tables with two real studies that disagree on GENE, and one simulated record."""
    data = tmp_path / "site"
    data.mkdir()
    pd.DataFrame([
        _evidence_row("E1", "REAL_A", GENE, "up", 0.01),
        _evidence_row("E2", "REAL_B", GENE, "down", 0.02),
        _evidence_row("E3", "REAL_A", "LOC900002", "hyper", 0.04, feature_type="methylation_region"),
        _evidence_row("E4", ALL_DEMO_STUDY_IDS[0], GENE, "up", 0.001, simulated=True),
    ]).to_parquet(data / "evidence.parquet", index=False)
    pd.DataFrame([
        {"feature_id_standardized": GENE, "feature_type": "gene", "phenotype": "thermal_tolerance",
         "simulated": False, "species_taxid": 29159, "tier": "high_priority_cross_study",
         "evidence_class": "resilience_associated", "contrast_type": "exposure_response",
         "is_predictive_baseline": False, "score": 0.7, "k_studies": 2, "studies": "REAL_A|REAL_B",
         "card_file": f"{GENE}_thermal.md"},
        {"feature_id_standardized": GENE, "feature_type": "gene", "phenotype": "larval_viability",
         "simulated": True, "species_taxid": 29159, "tier": "emerging", "evidence_class": "stress_response",
         "contrast_type": "exposure_response", "is_predictive_baseline": False, "score": 0.9,
         "k_studies": 1, "studies": ALL_DEMO_STUDY_IDS[0], "card_file": None},
    ]).to_parquet(data / "candidates.parquet", index=False)
    pd.DataFrame([{"card_file": f"{GENE}_thermal.md", "markdown": f"# {GENE} evidence card"}]).to_parquet(
        data / "cards.parquet", index=False)
    (data / "meta.json").write_text(json.dumps({"generated": "2026-10-01T00:00:00+00:00"}))
    return data


@pytest.fixture
def sources(site):
    return query.Sources.from_site_data(site)


# ---------------------------------------------------------------- a scripted Claude


def _usage(input_tokens=1000, output_tokens=200):
    return SimpleNamespace(input_tokens=input_tokens, output_tokens=output_tokens,
                           cache_creation_input_tokens=0, cache_read_input_tokens=0)


def text_reply(text, stop_reason="end_turn"):
    return {"texts": [text], "content": [SimpleNamespace(type="text", text=text)], "stop_reason": stop_reason}


def tool_reply(name, tool_input, tool_id="toolu_1", stop_reason="tool_use", preamble="Let me check. "):
    return {"texts": [preamble],
            "content": [SimpleNamespace(type="text", text=preamble),
                        SimpleNamespace(type="tool_use", id=tool_id, name=name, input=tool_input)],
            "stop_reason": stop_reason}


class FakeStream:
    def __init__(self, reply):
        self.reply = reply

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def __iter__(self):
        if isinstance(self.reply, Exception):
            raise self.reply
        for text in self.reply["texts"]:
            yield SimpleNamespace(type="text", text=text)

    def get_final_message(self):
        return SimpleNamespace(content=self.reply["content"], stop_reason=self.reply["stop_reason"],
                               usage=_usage(), model="claude-opus-5-5")


class FakeClient:
    """Answers each request with the next scripted reply, recording what it was sent."""

    def __init__(self, *replies):
        self.replies = list(replies)
        self.calls = []
        self.beta = SimpleNamespace(messages=SimpleNamespace(stream=self._stream))

    def _stream(self, **params):
        self.calls.append(copy.deepcopy(params))
        return FakeStream(self.replies.pop(0))


def _turn(client, sources, question="Is LOC900001 linked to thermal tolerance?", config=None, conv=None,
          spend=lambda usd: None):
    conv = conv or Conversation(id="c")
    events = list(run_turn(client, conv, question, sources, config or ChatConfig(), "SYSTEM", spend))
    return events, conv


# ---------------------------------------------------------------- tools


def test_tools_serve_real_studies_only_with_compact_records(sources):
    content, is_error = tools.run_tool("get_evidence", {"feature": GENE}, sources)
    out = json.loads(content)
    assert not is_error
    assert {r["study_id"] for r in out["records"]} == {"REAL_A", "REAL_B"}  # the simulated record is gone
    assert out["conflicting_directions"] is True
    assert set(out["records"][0]) == set(tools.CHAT_EVIDENCE_COLUMNS)
    assert "artifacts" not in out

    candidates = json.loads(tools.run_tool("list_candidates", {}, sources)[0])
    assert [c["phenotype"] for c in candidates["candidates"]] == ["thermal_tolerance"]
    assert "card_file" not in candidates["candidates"][0]

    card = json.loads(tools.run_tool("get_evidence_card", {"feature": GENE}, sources)[0])
    assert card["cards"][0]["markdown"] == f"# {GENE} evidence card"

    demo = json.loads(tools.run_tool("get_study", {"study_id": ALL_DEMO_STUDY_IDS[0]}, sources)[0])
    assert "simulated demo study" in demo["error"]
    missing = json.loads(tools.run_tool("get_study", {"study_id": "NOPE"}, sources)[0])
    assert not set(missing["known_study_ids"]) & set(ALL_DEMO_STUDY_IDS)


def test_replication_counts_real_studies_only(sources):
    out = json.loads(tools.run_tool("find_replicated_features", {"phenotype": "thermal_tolerance"}, sources)[0])
    (hit,) = out["features"]
    assert hit["feature_id_standardized"] == GENE
    assert hit["studies"] == ["REAL_A", "REAL_B"]  # the simulated study's record is not counted
    assert (hit["n_studies_up"], hit["n_studies_down"], hit["conflicting_directions"]) == (1, 1, True)
    one_study = json.loads(tools.run_tool("find_replicated_features", {"min_studies": 1}, sources)[0])
    assert "error" in one_study  # replication needs at least two studies


def test_a_shortened_result_says_how_much_was_left_out(sources):
    partial = json.loads(tools.run_tool("get_evidence", {"feature": GENE, "limit": 1}, sources)[0])
    assert partial["n_records"] == 2 and len(partial["records"]) == 1
    assert "first 1 of 2 records" in partial["partial_result"]
    whole = json.loads(tools.run_tool("get_evidence", {"feature": GENE}, sources)[0])
    assert "partial_result" not in whole


def test_tool_inputs_are_validated_before_anything_runs(sources):
    for name, args in [("get_evidence", {"limit": 5000}), ("get_evidence", {"drop_table": "x"}),
                       ("get_study", {}), ("list_vocabulary", {"name": "../secrets"}), ("no_such_tool", {})]:
        content, is_error = tools.run_tool(name, args, sources)
        assert is_error and "error" in json.loads(content)


def test_every_tool_streams_its_input_and_runs(sources):
    calls = {"aree_overview": {}, "list_vocabulary": {"name": "phenotype"}, "list_studies": {},
             "get_study": {"study_id": "NOPE"}, "search_features": {"query_text": "LOC9"},
             "get_evidence": {"phenotype": "thermal_tolerance"}, "list_candidates": {"tier": "emerging"},
             "find_replicated_features": {"phenotype": "thermal_tolerance"},
             "get_evidence_card": {"feature": GENE}}
    assert {t["name"] for t in tools.TOOLS} == set(calls)
    assert all(t["eager_input_streaming"] for t in tools.TOOLS)
    for name, args in calls.items():
        json.loads(tools.run_tool(name, args, sources)[0])


def test_oversized_results_are_cut_to_fit_and_say_so():
    big = {"n_records": 5000, "records": [{"evidence_id": f"E{i}", "pad": "x" * 200} for i in range(5000)]}
    text = tools._fit(big)
    assert len(text) <= tools.MAX_RESULT_CHARS
    out = json.loads(text)
    assert out["n_records"] == 5000 and len(out["records"]) < 5000
    assert "full query" in out["shortened_for_chat"]


# ---------------------------------------------------------------- one turn


def test_a_turn_runs_the_tool_and_answers_from_its_result(sources):
    client = FakeClient(tool_reply("get_evidence", {"feature": GENE}),
                        text_reply("Two real studies disagree (REAL_A up, REAL_B down)."))
    spent = []
    events, conv = _turn(client, sources, spend=spent.append)

    assert [e["type"] for e in events] == ["text", "tool_call", "tool_result", "text", "done"]
    assert events[1]["name"] == "get_evidence" and events[1]["input"] == {"feature": GENE}
    assert json.loads(events[2]["content"])["n_records"] == 2
    done = events[-1]
    assert done["tools"] == ["get_evidence"]
    assert len(spent) == 2 and done["cost_usd"] == pytest.approx(sum(spent)) and sum(spent) > 0

    first, second = client.calls
    assert first["system"][0]["cache_control"] == {"type": "ephemeral"}
    assert first["cache_control"] == {"type": "ephemeral"}
    assert first["fallbacks"] == "default" and first["betas"] == ["server-side-fallback-2026-07-01"]
    assert first["output_config"] == {"effort": "medium"} and "tool_choice" not in first
    # Claude's reply goes back unchanged, followed by the tool result it asked for.
    result = second["messages"][-1]["content"][0]
    assert result["tool_use_id"] == "toolu_1" and result["content"] == events[2]["content"]
    assert [m["role"] for m in conv.messages] == ["user", "assistant", "user", "assistant"]
    assert conv.turns == 1


def test_the_round_after_the_last_tool_round_must_answer(sources):
    client = FakeClient(tool_reply("aree_overview", {}), text_reply("Here is what AREE holds."))
    events, _ = _turn(client, sources, config=ChatConfig(max_tool_rounds=1))
    assert client.calls[1]["tool_choice"] == {"type": "none"}
    assert events[-1]["type"] == "done"


def test_a_truncated_tool_input_is_not_run(sources):
    client = FakeClient(tool_reply("get_evidence", {"feature": "LOC9"}, stop_reason="max_tokens"),
                        text_reply("Please narrow the question."))
    events, _ = _turn(client, sources)
    result = next(e for e in events if e["type"] == "tool_result")
    assert result["is_error"] and "cut off" in result["content"]


def test_a_refusal_is_reported_not_shown_as_an_answer(sources):
    client = FakeClient({"texts": [], "content": [], "stop_reason": "refusal"})
    events, _ = _turn(client, sources)
    assert events[0]["type"] == "error" and events[0]["code"] == "refusal"


def test_an_unparseable_tool_input_retries_the_round(sources):
    client = FakeClient(ValueError("bad JSON"), text_reply("Answer."))
    events, _ = _turn(client, sources)
    assert [e["type"] for e in events] == ["round_reset", "text", "done"]


def test_a_failed_turn_leaves_no_half_built_history(sources):
    conv = Conversation(id="c", messages=[{"role": "user", "content": "earlier"}])
    client = FakeClient(tool_reply("aree_overview", {}), RuntimeError("network"))
    with pytest.raises(RuntimeError):
        _turn(client, sources, conv=conv)
    assert conv.messages == [{"role": "user", "content": "earlier"}] and conv.turns == 0


def test_conversations_have_a_question_limit_and_expire():
    events = list(run_turn(FakeClient(), Conversation(id="c", turns=2), "q", None, ChatConfig(max_turns=2), "S"))
    assert events[0]["code"] == "conversation_limit"

    now = [0.0]
    store = Conversations(ttl_seconds=60, clock=lambda: now[0])
    first = store.get_or_create(None)
    assert store.get_or_create(first.id) is first
    now[0] = 120
    assert store.get_or_create(first.id) is not first


# ---------------------------------------------------------------- limits


def test_rate_limits_and_the_daily_budget():
    now = [0.0]
    limiter = RateLimiter(per_hour=2, per_day=3, clock=lambda: now[0])
    assert limiter.check("a") is None and limiter.check("a") is None
    assert "hour" in limiter.check("a")
    assert limiter.check("b") is None  # per client
    now[0] = 3601
    assert limiter.check("a") is None
    assert "day" in limiter.check("a")

    day = ["d1"]
    budget = DailyBudget(1.0, today=lambda: day[0])
    budget.add(0.6)
    assert not budget.exhausted()
    budget.add(0.6)
    assert budget.exhausted()
    day[0] = "d2"
    assert not budget.exhausted()

    usage = SimpleNamespace(input_tokens=1_000_000, output_tokens=1_000_000,
                            cache_creation_input_tokens=0, cache_read_input_tokens=1_000_000)
    assert cost_usd("claude-opus-5-5", usage) == pytest.approx(4.0 + 20.0 + 0.2)
    assert cost_usd("a-future-model", usage) >= cost_usd("claude-opus-5-5", usage)


# ---------------------------------------------------------------- web API


def _events(response) -> list[dict]:
    return [json.loads(line[len("data: "):]) for line in response.text.splitlines() if line.startswith("data: ")]


def _app(site, client, **overrides):
    settings = Settings(data_dir=site, **overrides)
    return TestClient(create_app(settings, client=client, verify=lambda secret, token, ip: token == "ok"))


def test_chat_streams_events_and_keeps_the_conversation(site):
    client = FakeClient(text_reply("First answer."), text_reply("Second answer."))
    app = _app(site, client)
    assert app.get("/api/health").json()["data_build"] == "2026-10-01T00:00:00+00:00"

    first = _events(app.post("/api/chat", json={"message": "What does AREE hold?"}))
    assert [e["type"] for e in first] == ["start", "text", "done"]
    conversation_id = first[0]["conversation_id"]
    second = _events(app.post("/api/chat", json={"message": "And for heat?", "conversation_id": conversation_id}))
    assert second[0]["conversation_id"] == conversation_id
    assert [m["role"] for m in client.calls[1]["messages"]] == ["user", "assistant", "user"]
    assert "published build of 2026-10-01" in client.calls[0]["system"][0]["text"]


def test_chat_refuses_what_its_limits_forbid(site):
    app = _app(site, FakeClient(*[text_reply("ok")] * 5), per_ip_hourly=2, max_question_chars=20)
    assert app.post("/api/chat", json={"message": "x" * 21}).status_code == 400
    assert app.post("/api/chat", json={"message": "   "}).status_code == 400
    assert app.post("/api/chat", json={"message": "one"}).status_code == 200
    assert app.post("/api/chat", json={"message": "two"}).status_code == 200
    refused = app.post("/api/chat", json={"message": "three"})
    assert refused.status_code == 429 and refused.json()["error"] == "rate_limit"

    broke = _app(site, FakeClient(text_reply("ok")), daily_budget_usd=0.0)
    assert broke.post("/api/chat", json={"message": "hi"}).status_code == 503

    guarded = _app(site, FakeClient(text_reply("ok")), turnstile_secret="s")
    assert guarded.post("/api/chat", json={"message": "hi", "turnstile_token": "forged"}).status_code == 403
    assert guarded.post("/api/chat", json={"message": "hi", "turnstile_token": "ok"}).status_code == 200


def test_only_the_site_and_local_pages_may_call_the_api(site):
    app = _app(site, FakeClient())
    allowed = app.options("/api/chat", headers={"Origin": "https://sr320.github.io",
                                                "Access-Control-Request-Method": "POST"})
    assert allowed.headers.get("access-control-allow-origin") == "https://sr320.github.io"
    other = app.options("/api/chat", headers={"Origin": "https://example.com",
                                              "Access-Control-Request-Method": "POST"})
    assert "access-control-allow-origin" not in other.headers
    local = _app(site, FakeClient(), allow_localhost=True).options(
        "/api/chat", headers={"Origin": "http://localhost:4200", "Access-Control-Request-Method": "POST"})
    assert local.headers.get("access-control-allow-origin") == "http://localhost:4200"


def test_settings_come_from_the_environment(monkeypatch, site):
    monkeypatch.setenv("AREE_DATA_DIR", str(site))
    monkeypatch.setenv("AREE_CHAT_MODEL", "claude-haiku-5-5")
    monkeypatch.setenv("AREE_CHAT_DAILY_BUDGET_USD", "2.5")
    monkeypatch.setenv("AREE_CHAT_ALLOWED_ORIGINS", "https://a.example, https://b.example")
    monkeypatch.setenv("AREE_CHAT_TRUST_PROXY", "true")
    settings = Settings.from_env()
    assert settings.data_dir == site and settings.daily_budget_usd == 2.5 and settings.trust_proxy
    assert settings.allowed_origins == ["https://a.example", "https://b.example"]
    assert settings.chat_config().fallbacks is False  # Haiku has no server-side fallback


def test_an_unexpected_failure_still_ends_the_stream_with_an_error(site):
    class Broken:
        beta = SimpleNamespace(messages=SimpleNamespace(stream=lambda **params: (_ for _ in ()).throw(
            TypeError("Could not resolve authentication method"))))

    events = _events(_app(site, Broken()).post("/api/chat", json={"message": "hi"}))
    assert [e["type"] for e in events] == ["start", "error"] and events[-1]["code"] == "internal"


def test_the_container_entry_point_serves_from_the_environment(monkeypatch, site):
    from chat.server import app_from_env

    monkeypatch.setenv("AREE_DATA_DIR", str(site))
    monkeypatch.setenv("ANTHROPIC_API_KEY", "test-key-not-used")
    health = TestClient(app_from_env()).get("/api/health").json()
    assert health["ok"] and health["data_build"] == "2026-10-01T00:00:00+00:00"
