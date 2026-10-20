"""Register validated studies into the flat registry/study_registry.csv index.

The CSV is a generated index, not a source of truth: the study YAMLs are.
It is gitignored so that study branches developed in parallel never conflict
on it; `rebuild_registry` regenerates it from every YAML in one pass.

Simulated demo studies (registry/demo_studies/) are indexed separately, in
registry/demo_study_registry.csv, so the study registry lists real studies only.
"""
from __future__ import annotations

import csv

from common import DEMO_REGISTRY_CSV, DEMO_STUDIES_DIR, STUDIES_DIR, STUDY_REGISTRY_CSV, load_yaml

from .schema_validate import validate_study_file

REGISTRY_FIELDS = [
    "study_id", "species", "species_identity_basis", "strain_or_population", "genome_assembly", "annotation_version",
    "assay_type", "analysis_mode", "platform", "raw_data_available", "processed_data_available",
    "data_status", "qc_status", "analysis_status", "n_comparisons", "simulated", "citation",
    "registered_by", "date_registered",
]


class DuplicateStudyError(RuntimeError):
    pass


def _registry_csv(simulated: bool):
    return DEMO_REGISTRY_CSV if simulated else STUDY_REGISTRY_CSV


def _read_registry(simulated: bool = False) -> list[dict]:
    path = _registry_csv(simulated)
    if not path.exists():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _write_registry(rows: list[dict], simulated: bool = False) -> None:
    path = _registry_csv(simulated)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=REGISTRY_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _study_to_row(study: dict) -> dict:
    return {
        "study_id": study["study_id"],
        "species": study.get("species", ""),
        # How the species was established, for studies in a hybridizing complex
        # (docs/species_identity.md); empty where the YAML records none.
        "species_identity_basis": (study.get("species_identity") or {}).get("basis", ""),
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

    A simulated study goes into the demo index (registry/demo_study_registry.csv)
    instead, so demo fixtures never appear in the study registry.

    Raises DuplicateStudyError if the study_id is already registered and
    allow_update is False. Raises ValueError if the study file fails schema
    validation.
    """
    result = validate_study_file(path)
    if not result.valid:
        raise ValueError(f"Study file failed validation: {'; '.join(result.errors)}")

    study = load_yaml(path)
    study_id = study["study_id"]
    simulated = bool(study.get("simulated", False))

    rows = _read_registry(simulated)
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

    _write_registry(rows, simulated)
    return row


def _index_dir(studies_dir) -> tuple[list[dict], dict[str, list[str]]]:
    rows, invalid = [], {}
    for path in sorted(studies_dir.glob("*.yaml")):
        if path.name.startswith("_"):
            continue
        result = validate_study_file(path)
        if not result.valid:
            invalid[str(path)] = result.errors
            continue
        rows.append(_study_to_row(load_yaml(path)))
    return rows, invalid


def rebuild_registry(studies_dir=None, demo_studies_dir=None) -> tuple[list[dict], dict[str, list[str]]]:
    """Regenerate both registry CSVs from every study YAML, replacing what is there.

    Real studies (registry/studies/) go to the study registry and simulated demo
    studies (registry/demo_studies/) to the demo index; the returned rows are the
    real ones. Files starting with `_` (templates) are skipped. A study that fails
    validation is left out of the index and reported in the returned
    {path: errors} map rather than aborting the rebuild; the caller decides
    whether that is fatal.
    """
    rows, invalid = _index_dir(studies_dir or STUDIES_DIR)
    demo_rows, demo_invalid = _index_dir(demo_studies_dir or DEMO_STUDIES_DIR)
    _write_registry(rows)
    _write_registry(demo_rows, simulated=True)
    return rows, {**invalid, **demo_invalid}


def list_studies(include_simulated: bool = False) -> list[dict]:
    """Rows of the study registry; with include_simulated, the demo index rows follow."""
    rows = _read_registry()
    if include_simulated:
        rows += _read_registry(simulated=True)
    return rows
