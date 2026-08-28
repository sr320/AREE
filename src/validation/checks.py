"""Small reusable validation checks, exercised directly by tests/ and by
src/intake for CLI-facing validation."""
from __future__ import annotations

from intake.schema_validate import ValidationResult, validate_study_file

REQUIRED_PROVENANCE_FIELDS = ["registered_by", "date_registered"]


def validate_study_schema(path) -> ValidationResult:
    return validate_study_file(path)


def check_duplicate_study_id(existing_ids, new_id: str) -> bool:
    """Returns True if new_id would be a duplicate of an existing registered study_id."""
    return new_id in set(existing_ids)


def check_required_provenance_fields(study: dict) -> list:
    """Returns a list of missing required provenance field names (empty if all present)."""
    provenance = study.get("provenance") or {}
    return [f for f in REQUIRED_PROVENANCE_FIELDS if not provenance.get(f)]
