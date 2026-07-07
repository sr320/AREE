"""Register validated studies into the flat registry/study_registry.csv index."""
from __future__ import annotations

import csv

from common import STUDY_REGISTRY_CSV, load_yaml

from .schema_validate import validate_study_file

REGISTRY_FIELDS = [
    "study_id", "species", "strain_or_population", "genome_assembly", "annotation_version",
    "assay_type", "analysis_mode", "platform", "raw_data_available", "processed_data_available",
    "data_status", "qc_status", "analysis_status", "n_comparisons", "simulated", "citation",
    "registered_by", "date_registered",
]


class DuplicateStudyError(RuntimeError):
    pass


def _read_registry() -> list[dict]:
    if not STUDY_REGISTRY_CSV.exists():
        return []
    with open(STUDY_REGISTRY_CSV, newline="") as fh:
        return list(csv.DictReader(fh))


def _write_registry(rows: list[dict]) -> None:
    STUDY_REGISTRY_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(STUDY_REGISTRY_CSV, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=REGISTRY_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _study_to_row(study: dict) -> dict:
    return {
        "study_id": study["study_id"],
        "species": study.get("species", ""),
        "strain_or_population": study.get("strain_or_population") or "",
        "genome_assembly": study.get("genome_assembly", ""),
        "annotation_version": study.get("annotation_version") or "",
        "assay_type": "|".join(study.get("assay_type", [])),
        "analysis_mode": study.get("analysis_mode", ""),
        "platform": study.get("platform") or "",
        "raw_data_available": study["data_availability"]["raw_data_available"],
        "processed_data_available": study["data_availability"]["processed_data_available"],
        "data_status": study["data_availability"]["status"],
        "qc_status": study.get("qc_status", "not_started"),
        "analysis_status": study.get("analysis_status", "not_started"),
        "n_comparisons": len(study.get("comparisons", [])),
        "simulated": study.get("simulated", False),
        "citation": study.get("citation", ""),
        "registered_by": study["provenance"]["registered_by"],
        "date_registered": study["provenance"]["date_registered"],
    }


def register_study(path, allow_update: bool = False) -> dict:
    """Validate and register a study YAML file into the study registry CSV.

    Raises DuplicateStudyError if the study_id is already registered and
    allow_update is False. Raises ValueError if the study file fails schema
    validation.
    """
    result = validate_study_file(path)
    if not result.valid:
        raise ValueError(f"Study file failed validation: {'; '.join(result.errors)}")

    study = load_yaml(path)
    study_id = study["study_id"]

    rows = _read_registry()
    existing_idx = next((i for i, r in enumerate(rows) if r["study_id"] == study_id), None)

    if existing_idx is not None and not allow_update:
        raise DuplicateStudyError(
            f"study_id {study_id!r} is already registered. Pass allow_update=True "
            "(or --update on the CLI) to overwrite the existing entry."
        )

    row = _study_to_row(study)
    if existing_idx is not None:
        rows[existing_idx] = row
    else:
        rows.append(row)

    _write_registry(rows)
    return row


def list_studies() -> list[dict]:
    return _read_registry()
