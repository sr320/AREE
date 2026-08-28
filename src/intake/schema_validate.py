"""Validate a study registration file against schemas/study.schema.json and
the controlled vocabularies referenced by its fields.

Errors are hard failures (schema violations, unknown vocabulary terms).
Warnings are soft signals (e.g. resilience_classification does not match the
ontology's default resilience_relevance for the chosen phenotype) that a
curator may have deliberately overridden with a documented rationale.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import jsonschema

from common import SCHEMAS_DIR, load_json, load_vocab, load_yaml


@dataclass
class ValidationResult:
    path: str
    valid: bool
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    def __bool__(self):
        return self.valid


def _load_study_schema() -> dict:
    return load_json(SCHEMAS_DIR / "study.schema.json")


def validate_study_file(path) -> ValidationResult:
    path = Path(path)
    errors: list[str] = []
    warnings: list[str] = []

    try:
        study = load_yaml(path)
    except Exception as exc:  # malformed YAML
        return ValidationResult(str(path), False, errors=[f"Could not parse YAML: {exc}"])

    schema = _load_study_schema()
    validator = jsonschema.Draft7Validator(schema)
    for err in sorted(validator.iter_errors(study), key=lambda e: e.path):
        loc = "/".join(str(p) for p in err.path) or "<root>"
        errors.append(f"{loc}: {err.message}")

    if errors:
        return ValidationResult(str(path), False, errors=errors, warnings=warnings)

    stem = path.stem
    if stem not in ("_TEMPLATE",) and study.get("study_id") != stem:
        errors.append(
            f"study_id {study.get('study_id')!r} does not match filename stem {stem!r}"
        )

    phenotype_vocab = load_vocab("phenotype_ontology")
    stressor_ids = set(load_vocab("stressor_ontology").keys())
    tissue_ids = set(load_vocab("tissue_types").keys())
    life_stage_ids = set(load_vocab("life_stages").keys())
    assay_ids = set(load_vocab("assay_types").keys())
    quality_flag_ids = set(load_vocab("quality_flags").keys())

    for assay in study.get("assay_type", []):
        if assay not in assay_ids:
            errors.append(f"assay_type {assay!r} is not in assay_types vocabulary")

    for flag in study.get("quality_flags", []) or []:
        if flag not in quality_flag_ids:
            errors.append(f"quality_flags entry {flag!r} is not in quality_flags vocabulary")

    seen_comparison_ids = set()
    for comp in study.get("comparisons", []):
        cid = comp.get("comparison_id")
        if cid in seen_comparison_ids:
            errors.append(f"duplicate comparison_id {cid!r} within study")
        seen_comparison_ids.add(cid)

        if comp.get("tissue") not in tissue_ids:
            errors.append(f"comparison {cid!r}: tissue {comp.get('tissue')!r} not in tissue_types vocabulary")
        if comp.get("life_stage") not in life_stage_ids:
            errors.append(f"comparison {cid!r}: life_stage {comp.get('life_stage')!r} not in life_stages vocabulary")
        if comp.get("stressor_standardized") not in stressor_ids:
            errors.append(
                f"comparison {cid!r}: stressor_standardized {comp.get('stressor_standardized')!r} "
                "not in stressor_ontology vocabulary"
            )
        phenotype_id = comp.get("phenotype")
        if phenotype_id not in phenotype_vocab:
            errors.append(f"comparison {cid!r}: phenotype {phenotype_id!r} not in phenotype_ontology vocabulary")
        else:
            expected_relevance = phenotype_vocab[phenotype_id]["resilience_relevance"]
            declared = comp.get("resilience_classification")
            if declared != expected_relevance:
                warnings.append(
                    f"comparison {cid!r}: resilience_classification={declared!r} differs from "
                    f"phenotype_ontology default ({expected_relevance!r}) for phenotype "
                    f"{phenotype_id!r}. Confirm this override is intentional and documented."
                )

    valid = len(errors) == 0
    return ValidationResult(str(path), valid, errors=errors, warnings=warnings)
