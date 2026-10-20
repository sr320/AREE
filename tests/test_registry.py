import shutil

import pytest

from common import DEMO_STUDIES_DIR, STUDIES_DIR, study_files, study_ids, study_yaml_path
from intake.registry import DuplicateStudyError, list_studies, register_study
from intake.schema_validate import validate_study_file
from validation.checks import check_duplicate_study_id, check_required_provenance_fields

REAL_STUDY = "registry/studies/HESSER2024_VCOR.yaml"
DEMO_STUDY = "registry/demo_studies/GIGAS_HEAT01.yaml"


def test_register_study_adds_to_registry(isolated_registry):
    register_study(REAL_STUDY)
    rows = list_studies()
    assert len(rows) == 1
    assert rows[0]["study_id"] == "HESSER2024_VCOR"


def test_register_study_twice_raises_duplicate_error(isolated_registry):
    register_study(REAL_STUDY)
    with pytest.raises(DuplicateStudyError):
        register_study(REAL_STUDY)


def test_register_study_allows_explicit_update(isolated_registry):
    register_study(REAL_STUDY)
    row = register_study(REAL_STUDY, allow_update=True)
    assert row["study_id"] == "HESSER2024_VCOR"
    assert len(list_studies()) == 1


def test_simulated_study_is_indexed_apart_from_the_study_registry(isolated_registry):
    register_study(REAL_STUDY)
    register_study(DEMO_STUDY)
    assert [r["study_id"] for r in list_studies()] == ["HESSER2024_VCOR"]
    assert [r["study_id"] for r in list_studies(include_simulated=True)] == ["HESSER2024_VCOR", "GIGAS_HEAT01"]
    with pytest.raises(DuplicateStudyError):
        register_study(DEMO_STUDY)


def test_real_and_demo_studies_live_in_their_own_directories():
    for path in study_files(include_demo=False):
        assert path.parent == STUDIES_DIR
    real = set(study_ids(simulated=False))
    demo = set(study_ids(simulated=True))
    assert demo and real and not real & demo
    assert {p.stem for p in DEMO_STUDIES_DIR.glob("*.yaml")} == demo
    assert study_yaml_path("GIGAS_HEAT01") == DEMO_STUDIES_DIR / "GIGAS_HEAT01.yaml"
    assert study_yaml_path("HESSER2024_VCOR") == STUDIES_DIR / "HESSER2024_VCOR.yaml"


def test_study_in_the_wrong_directory_fails_validation(tmp_path, monkeypatch):
    import intake.schema_validate as sv

    real_dir, demo_dir = tmp_path / "studies", tmp_path / "demo_studies"
    real_dir.mkdir()
    demo_dir.mkdir()
    monkeypatch.setattr(sv, "STUDIES_DIR", real_dir)
    monkeypatch.setattr(sv, "DEMO_STUDIES_DIR", demo_dir)

    misplaced_demo = shutil.copy(DEMO_STUDY, real_dir)
    misplaced_real = shutil.copy(REAL_STUDY, demo_dir)
    assert any("registry/demo_studies/" in e for e in validate_study_file(misplaced_demo).errors)
    assert any("registry/studies/" in e for e in validate_study_file(misplaced_real).errors)


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
