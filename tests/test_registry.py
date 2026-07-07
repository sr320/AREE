import pytest

from intake.registry import DuplicateStudyError, list_studies, register_study
from validation.checks import check_duplicate_study_id, check_required_provenance_fields


def test_register_study_adds_to_registry(isolated_registry):
    register_study("registry/studies/GIGAS_HEAT01.yaml")
    rows = list_studies()
    assert len(rows) == 1
    assert rows[0]["study_id"] == "GIGAS_HEAT01"


def test_register_study_twice_raises_duplicate_error(isolated_registry):
    register_study("registry/studies/GIGAS_HEAT01.yaml")
    with pytest.raises(DuplicateStudyError):
        register_study("registry/studies/GIGAS_HEAT01.yaml")


def test_register_study_allows_explicit_update(isolated_registry):
    register_study("registry/studies/GIGAS_HEAT01.yaml")
    row = register_study("registry/studies/GIGAS_HEAT01.yaml", allow_update=True)
    assert row["study_id"] == "GIGAS_HEAT01"
    assert len(list_studies()) == 1


def test_check_duplicate_study_id_helper():
    existing = ["GIGAS_HEAT01", "GIGAS_OA02"]
    assert check_duplicate_study_id(existing, "GIGAS_HEAT01") is True
    assert check_duplicate_study_id(existing, "GIGAS_NEW99") is False


def test_check_required_provenance_fields_flags_missing():
    complete = {"provenance": {"registered_by": "curator", "date_registered": "2026-01-01"}}
    incomplete = {"provenance": {"registered_by": "curator", "date_registered": ""}}
    missing_block = {}

    assert check_required_provenance_fields(complete) == []
    assert "date_registered" in check_required_provenance_fields(incomplete)
    assert set(check_required_provenance_fields(missing_block)) == {"registered_by", "date_registered"}
