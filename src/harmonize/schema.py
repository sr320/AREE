"""Canonical evidence-record column definition and shared row-building helpers.

Mirrors schemas/evidence.schema.json. Keeping the column list here (rather than
re-deriving it from the JSON schema at runtime) keeps DataFrame column order
stable and makes it obvious in code review when the two definitions diverge.
"""
from __future__ import annotations

import re

from common import resolve_species, sha256_file

from .identifiers import active_crosswalk_path, crosswalk_annotation_release

EVIDENCE_COLUMNS = [
    # `simulated` separates fabricated demo evidence from evidence derived from a
    # real published study. It is part of the schema rather than a report-time
    # annotation because the two must never be pooled, and a reader of the
    # evidence table must be able to tell them apart without consulting the registry.
    "evidence_id", "study_id", "comparison_id", "simulated",
    "feature_id_original", "feature_id_standardized", "feature_type",
    "orthogroup_id", "species", "species_taxid", "genome_assembly", "annotation_version",
    "identifier_annotation_release",
    "annotation_context", "molecular_direction", "effect_size", "effect_size_type",
    "standard_error", "ci_lower", "ci_upper", "p_value", "adjusted_p_value",
    "sample_size", "tissue", "life_stage", "stressor", "phenotype",
    "phenotype_direction", "analysis_method", "mapping_confidence",
    "quality_flags", "source_file", "workflow_version", "date_generated", "generated_by",
]

_SLUG_RE = re.compile(r"[^A-Za-z0-9_.-]+")


def _slug(value: str) -> str:
    return _SLUG_RE.sub("_", str(value))


def make_evidence_id(study_id: str, comparison_id: str, feature_id_original: str, analysis_method: str) -> str:
    return "__".join(_slug(v) for v in (study_id, comparison_id, feature_id_original, analysis_method))


def compute_quality_flags(study: dict, comparison: dict, mapping_confidence: str) -> list[str]:
    flags = set(study.get("quality_flags") or [])
    if comparison.get("biological_replicates", 0) < 3:
        flags.add("low_replication")
    if mapping_confidence in {"inferred", "one_to_many_ortholog", "many_to_one_ortholog", "unresolved"}:
        flags.add("identifier_mapping_uncertain")
    return sorted(flags)


def source_file_ref(path) -> str:
    return f"{path}#sha256:{sha256_file(path)}"


def molecular_direction_from_effect(effect_size: float | None, up_label="up", down_label="down") -> str:
    if effect_size is None:
        return "ambiguous"
    if effect_size > 0:
        return up_label
    if effect_size < 0:
        return down_label
    return "no_change"


def study_reference_fields(study: dict) -> dict:
    """Species/assembly/annotation columns shared by every evidence record.

    Centralized so the four assay harmonizers cannot drift apart on how a study's
    reference context is recorded.

    `species` preserves what the study reported; `species_taxid` is the canonical
    identity that meta-analysis groups on, so Crassostrea gigas and Magallana
    gigas do not become two species. `identifier_annotation_release` records the
    annotation the crosswalk resolved identifiers against, which is often a
    *different* assembly from `genome_assembly` — NCBI Gene IDs persist across
    assemblies, but a reader must be able to see that the crossing happened.
    """
    species = resolve_species(study.get("species"))
    return {
        "species": study["species"],
        "species_taxid": species.ncbi_taxid if species else None,
        "genome_assembly": study["genome_assembly"],
        "annotation_version": study.get("annotation_version"),
        "identifier_annotation_release": crosswalk_annotation_release(active_crosswalk_path()),
    }
