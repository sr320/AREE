"""AREE's query layer as Claude tools, for the public chat.

The same read-only queries as the MCP server (`query.mcp_server`), with
three differences suited to a public, pay-per-token chat:

* answers come from real studies only: simulated demo data never reaches a visitor;
* result sizes are capped, and evidence records carry a compact column set
  (an evidence_id leads back to the full provenance);
* every tool input is validated against its schema before anything runs, since
  the input is model output and the chat is open to anyone.
"""
from __future__ import annotations

import json
import logging

import jsonschema

import query

log = logging.getLogger(__name__)

# About 10k tokens: enough for a page of records or one evidence card.
MAX_RESULT_CHARS = 40_000

CHAT_EVIDENCE_COLUMNS = [
    "evidence_id", "study_id", "comparison_id", "contrast_type", "species",
    "feature_id_original", "feature_id_standardized", "feature_type", "orthogroup_id",
    "mapping_confidence", "molecular_direction", "effect_size", "effect_size_type",
    "standard_error", "ci_lower", "ci_upper", "adjusted_p_value", "sample_size",
    "tissue", "life_stage", "stressor", "phenotype", "phenotype_direction",
    "resilience_classification", "quality_flags",
]

_TEXT = {"type": "string", "minLength": 1, "maxLength": 200}


def _limit(default: int, maximum: int) -> dict:
    return {"type": "integer", "minimum": 1, "maximum": maximum,
            "description": f"How many rows to return (default {default}, at most {maximum})."}


def _schema(properties: dict, required: tuple[str, ...] = ()) -> dict:
    return {"type": "object", "properties": properties, "required": list(required),
            "additionalProperties": False}


TOOLS = [
    {
        "name": "aree_overview",
        "description": "Start here. Counts of studies, evidence records and candidates, and the "
                       "phenotypes, stressors, feature types and species that have evidence.",
        "input_schema": _schema({}),
    },
    {
        "name": "list_vocabulary",
        "description": "Controlled-vocabulary term ids and labels, to turn the visitor's words into "
                       "filter values (e.g. 'heat stress' -> stressor 'temperature').",
        "input_schema": _schema({"name": {"type": "string", "enum": sorted(query.VOCABULARIES)}}, ("name",)),
    },
    {
        "name": "list_studies",
        "description": "Registered studies, optionally filtered. stressor, phenotype and assay_type are "
                       "controlled ids (e.g. temperature, thermal_tolerance, rnaseq); species matches any "
                       "part of the name.",
        "input_schema": _schema({"species": _TEXT, "stressor": _TEXT, "phenotype": _TEXT,
                                 "assay_type": _TEXT}),
    },
    {
        "name": "get_study",
        "description": "One study's full registration record: accessions, citation, each comparison's "
                       "tissue, life stage, original and standardized stressor, phenotype, sample size, "
                       "limitations, QC, and how many evidence records it contributed.",
        "input_schema": _schema({"study_id": _TEXT}, ("study_id",)),
    },
    {
        "name": "search_features",
        "description": "Find genes, proteins, methylation regions or metabolites whose identifier "
                       "(standardized, original, or orthogroup) contains query_text, e.g. 'LOC105333935' or "
                       "'OG0001234'. One row per feature with its studies, phenotypes and directions.",
        "input_schema": _schema({"query_text": _TEXT, "limit": _limit(10, 25)}, ("query_text",)),
    },
    {
        "name": "get_evidence",
        "description": "Harmonized evidence records (effect size, SE/CI, adjusted p, direction, context, "
                       "mapping confidence, quality flags), most significant first, with the total count, "
                       "the number of studies, direction counts, and whether directions conflict. Filter "
                       "values are controlled ids; feature matches an identifier exactly unless "
                       "exact_feature is false.",
        "input_schema": _schema({
            "feature": _TEXT, "study_id": _TEXT, "phenotype": _TEXT, "stressor": _TEXT,
            "species": _TEXT, "feature_type": _TEXT, "tissue": _TEXT,
            "max_adjusted_p": {"type": "number", "minimum": 0, "maximum": 1},
            "exact_feature": {"type": "boolean"},
            "limit": _limit(20, 50),
        }),
    },
    {
        "name": "find_replicated_features",
        "description": "Features significant (adjusted p <= max_adjusted_p, default 0.05) in at least "
                       "min_studies studies (default 2). Use it for questions like 'which genes are "
                       "significant for thermal_tolerance in more than one study, and do any disagree?'. "
                       "Per feature: the studies, how many went up vs down, whether directions conflict, "
                       "and the smallest adjusted p; plus how many matching features conflict overall.",
        "input_schema": _schema({
            "phenotype": _TEXT, "stressor": _TEXT, "feature_type": _TEXT, "species": _TEXT,
            "tissue": _TEXT, "min_studies": {"type": "integer", "minimum": 2, "maximum": 50},
            "max_adjusted_p": {"type": "number", "minimum": 0, "maximum": 1},
            "limit": _limit(15, 25),
        }),
    },
    {
        "name": "list_candidates",
        "description": "Ranked candidate biomarkers with transparent score components (k_studies, "
                       "direction consistency, I-squared, pooled effect, adjusted p). tier: "
                       "high_priority_cross_study, multi_omics_convergence, metabolite_linked_convergence, "
                       "emerging. evidence_class: resilience_associated, stress_response, "
                       "disease_associated, exposure_only. contrast_type: baseline_phenotype_difference "
                       "(resilient vs susceptible animals under unstressed conditions, the design for a "
                       "predictive marker), phenotype_difference_under_exposure, exposure_response, "
                       "phenotype_by_exposure_interaction.",
        "input_schema": _schema({
            "phenotype": _TEXT, "tier": _TEXT, "evidence_class": _TEXT, "feature": _TEXT,
            "species_taxid": {"type": "integer", "minimum": 1}, "feature_type": _TEXT,
            "contrast_type": _TEXT, "limit": _limit(15, 25),
        }),
    },
    {
        "name": "get_evidence_card",
        "description": "The evidence card for a candidate: supporting studies, effect summary, contexts, "
                       "limitations, what the evidence does and does not show, and the recommended next "
                       "validation step.",
        "input_schema": _schema({"feature": _TEXT, "phenotype": _TEXT}, ("feature",)),
    },
]

for _tool in TOOLS:
    # Inputs stream as they are generated; run_tool validates them, as eager
    # streaming leaves validation to the caller.
    _tool["eager_input_streaming"] = True

_SCHEMAS = {tool["name"]: tool["input_schema"] for tool in TOOLS}


def _evidence(sources, **args):
    out = query.get_evidence(**{"limit": 20, **args}, include_simulated=False, sources=sources)
    out["records"] = [{k: r.get(k) for k in CHAT_EVIDENCE_COLUMNS} for r in out.get("records", [])]
    return out


def _candidates(sources, **args):
    out = query.list_candidates(**{"limit": 15, **args}, include_simulated=False, sources=sources)
    for rec in out.get("candidates", []):
        rec.pop("card_file", None)
    return out


def _study(sources, study_id):
    out = query.get_study(study_id, sources)
    if "study" in out and out["study"].get("simulated"):
        return {"error": f"{study_id} is a simulated demo study, not part of AREE's evidence."}
    if "known_study_ids" in out:
        real = {s["study_id"] for s in query.list_studies(include_simulated=False, sources=sources)["studies"]}
        out["known_study_ids"] = [s for s in out["known_study_ids"] if s in real]
    return out


_HANDLERS = {
    "aree_overview": lambda sources: query.overview(sources),
    "list_vocabulary": lambda sources, name: query.list_vocabulary(name),
    "list_studies": lambda sources, **args: query.list_studies(**args, include_simulated=False, sources=sources),
    "get_study": _study,
    "search_features": lambda sources, query_text, limit=10: query.search_features(
        query_text, limit, include_simulated=False, sources=sources),
    "get_evidence": _evidence,
    "find_replicated_features": lambda sources, **args: query.find_replicated_features(
        **{"limit": 15, **args}, include_simulated=False, sources=sources),
    "list_candidates": _candidates,
    "get_evidence_card": lambda sources, feature, phenotype=None: query.get_evidence_card(
        feature, phenotype, include_simulated=False, sources=sources),
}


def _fit(result: dict) -> str:
    """JSON for the model, cut to MAX_RESULT_CHARS by dropping rows from the end of the
    longest list (or shortening long text), and saying so."""
    result = {k: v for k, v in result.items() if k != "artifacts"}
    text = json.dumps(result, default=str)
    while len(text) > MAX_RESULT_CHARS:
        lists = [k for k, v in result.items() if isinstance(v, list) and v]
        if lists:
            key = max(lists, key=lambda k: len(json.dumps(result[k], default=str)))
            if len(result[key]) > 1:
                result[key] = result[key][: len(result[key]) // 2]
                result["shortened_for_chat"] = (f"Only the first {len(result[key])} {key} are shown; "
                                                "counts above are for the full query.")
                text = json.dumps(result, default=str)
                continue
        return text[:MAX_RESULT_CHARS] + ' ... [cut to fit; see the Explore page for the full result]'
    return text


def run_tool(name: str, args, sources: query.Sources) -> tuple[str, bool]:
    """Run one tool call; returns (JSON content for the model, whether it is an error)."""
    schema = _SCHEMAS.get(name)
    if schema is None:
        return json.dumps({"error": f"Unknown tool {name!r}."}), True
    try:
        jsonschema.validate(args, schema)
    except jsonschema.ValidationError as exc:
        return json.dumps({"error": f"Invalid input for {name}: {exc.message}"}), True
    try:
        result = _HANDLERS[name](sources, **args)
    except Exception:
        log.exception("tool %s failed with input %r", name, args)
        return json.dumps({"error": f"{name} failed on the server; try a narrower question."}), True
    return _fit(result), "error" in result
