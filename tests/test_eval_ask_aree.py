"""The Ask AREE eval's own machinery: its case file and its programmatic grader.

No API calls: the judge is replaced by a stub. These check that a known-good answer
passes, a known-bad one fails, and the case file stays consistent with the runner.
"""
from __future__ import annotations

import importlib.util

import yaml
from conftest import ALL_DEMO_STUDY_IDS

from common import REPO_ROOT

spec = importlib.util.spec_from_file_location("eval_ask_aree", REPO_ROOT / "scripts" / "eval_ask_aree.py")
ev = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ev)

REAL = {"LI2018_HEAT_TIDAL", "LUTIER2022_OA_TIPPING"}
SIMULATED = set(ALL_DEMO_STUDY_IDS)


def test_every_case_is_well_formed_and_uses_known_facts():
    cases = yaml.safe_load(ev.CASES_PATH.read_text())
    ids = [c["id"] for c in cases]
    assert len(ids) == len(set(ids))
    for case in cases:
        assert case["tags"] and case["turns"] and all(isinstance(t, str) and t.strip() for t in case["turns"])
        assert case.get("facts") or case.get("judge"), f"{case['id']} has nothing to grade it by"
        assert set(case.get("facts", [])) <= set(ev.FACTS), f"{case['id']} names an unknown fact"
    assert set(ev.HARNESS_PATHS) >= {"scripts/eval_ask_aree.py", "evals/ask_aree/cases.yaml"}
    assert all((REPO_ROOT / p).exists() for p in ev.HARNESS_PATHS)


def test_numbers_must_be_stated_exactly_unless_large_and_rounded():
    assert ev.check_fact("2,919 genes are significant", ("numbers", [2919]))[0]
    assert ev.check_fact("about 1.24 million records", ("numbers", [1236316]))[0]
    assert not ev.check_fact("about 2.9 thousand genes", ("numbers", [2919]))[0]
    assert not ev.check_fact("291 genes", ("numbers", [2919]))[0]
    assert not ev.check_fact("", ("numbers", [35]))[0]
    assert ev.check_fact("either 30 or 31", ("numbers", [31, 99]))[0]  # any acceptable value


def test_mentions_match_ids_and_labels_in_any_form():
    phenotypes = ("all", ["thermal_tolerance", "disease_resistance"])
    assert ev.check_fact("Thermal tolerance and disease resistance.", phenotypes)[0]
    assert ev.check_fact("`thermal_tolerance`, `disease_resistance`", phenotypes)[0]
    ok, why = ev.check_fact("Only thermal tolerance.", phenotypes)
    assert not ok and "missing disease_resistance" in why
    assert ev.check_fact("the Manila clam", ("any", ["Ruditapes philippinarum", "Manila clam"]))[0]
    assert not ev.check_fact("I don't know.", ("any", ["Ruditapes philippinarum", "Manila clam"]))[0]


def test_grounding_flags_invented_and_simulated_study_ids_only():
    assert ev.grounded("See LI2018_HEAT_TIDAL and LUTIER2022_OA_TIPPING.", REAL, SIMULATED)[0]
    ok, why = ev.grounded("Supported by SMITH2021_HEAT_GILL.", REAL, SIMULATED)
    assert not ok and "unknown" in why
    demo = ALL_DEMO_STUDY_IDS[0]
    ok, why = ev.grounded(f"See {demo}.", REAL, SIMULATED)
    assert not ok and "simulated" in why
    # Evidence ids, assembly accessions, orthogroups and plain text are not study ids.
    assert ev.grounded("Record LI2018_HEAT_TIDAL__heat_35c_24h_vs_0h__109620656__DESeq2 on GCF_963853765.1 "
                       "in OG0000560; HIGH_PRIORITY tier.", REAL, SIMULATED)[0]


def test_a_good_answer_passes_and_a_bad_one_fails_end_to_end():
    case = {"id": "c", "turns": ["Which studies measured acidification tolerance?"], "tags": ["studies"],
            "facts": ["acid"], "judge": "Lists them."}
    facts = {"acid": ("all", ["LUTIER2022_OA_TIPPING"])}

    def judge(case, trace):
        verdict = "LUTIER2022" in trace[-1]["content"]
        return ({m: float(verdict) for m in ev.JUDGE_METRICS}, dict.fromkeys(ev.JUDGE_METRICS, "stub"),
                {"input_tokens": 1, "output_tokens": 1})

    good = "LUTIER2022_OA_TIPPING measured it."
    scores, why, _ = ev.grade(case, good, [{"role": "assistant", "content": good}], facts, (REAL, SIMULATED), judge)
    assert scores["pass"] == 1.0 and set(scores) == {m["id"] for m in ev.METRICS}

    for bad in ["", "I don't know.", "SMITH2021_OA_GILL measured it."]:
        scores, _, _ = ev.grade(case, bad, [{"role": "assistant", "content": bad}], facts, (REAL, SIMULATED), judge)
        assert scores["pass"] == 0.0 and scores["facts"] == 0.0


def test_the_judge_sees_tool_results_whole_as_the_chat_model_did():
    from chat.tools import MAX_RESULT_CHARS

    # A shorter limit made the judge call claims from the unseen part unsupported.
    assert ev.MAX_JUDGE_RESULT_CHARS >= MAX_RESULT_CHARS
    full = '{"features": [' + ", ".join(f'"G{i}"' for i in range(3000)) + "]}"
    assert len(full) < MAX_RESULT_CHARS
    text = ev.judge_input({"id": "c", "turns": ["q"]}, [{"role": "user", "content": "q"},
                                                       {"role": "tool_result", "content": full}])
    assert full in text


def test_the_judge_knows_what_the_chat_was_told():
    # Without the chat's instructions and tool descriptions, the judge called the build
    # date, the Explore link and contrast-type meanings "invented".
    context = ev.judge_context([{"role": "system", "content": "The data are AREE's published build of X."}])
    assert "published build of X" in context and "find_replicated_features" in context
    assert "baseline_phenotype_difference" in context
    assert "consistent before conflicting" in ev.JUDGE_SYSTEM


def test_the_judge_sees_every_tool_result_and_earlier_answers_as_data():
    case = {"id": "c", "turns": ["First?", "Then?"], "judge": "Be right."}
    trace = [{"role": "system", "content": "SYSTEM"}, {"role": "user", "content": "First?"},
             {"role": "tool_call", "name": "list_studies", "content": "{}"},
             {"role": "tool_result", "content": "x" * (ev.MAX_JUDGE_RESULT_CHARS + 50)},
             {"role": "assistant", "content": "Earlier answer."}, {"role": "user", "content": "Then?"},
             {"role": "assistant", "content": "Final answer."}]
    text = ev.judge_input(case, trace)
    assert "CALL list_studies" in text and "[result shortened for grading]" in text
    assert "<earlier_answers>\nEarlier answer." in text and "<answer>\nFinal answer.\n</answer>" in text
    assert "Be right." in text and "SYSTEM" not in text
