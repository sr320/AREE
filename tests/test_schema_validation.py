import glob
import textwrap

import yaml

from intake.schema_validate import validate_study_file


def test_all_demo_studies_are_valid():
    paths = sorted(glob.glob("registry/studies/GIGAS_*.yaml"))
    assert len(paths) == 6
    for path in paths:
        result = validate_study_file(path)
        assert result.valid, f"{path} failed validation: {result.errors}"


def test_malformed_yaml_is_rejected(tmp_path):
    bad_file = tmp_path / "BAD.yaml"
    bad_file.write_text("study_id: [unclosed\n  this: is not valid yaml")
    result = validate_study_file(bad_file)
    assert not result.valid
    assert result.errors


def test_missing_required_field_is_rejected(tmp_path):
    minimal = tmp_path / "MISSING_FIELDS.yaml"
    minimal.write_text(textwrap.dedent("""
        study_id: MISSING_FIELDS
        species: "Crassostrea gigas"
    """))
    result = validate_study_file(minimal)
    assert not result.valid
    assert result.errors


def test_unknown_vocabulary_term_is_rejected(tmp_path):
    with open("registry/studies/GIGAS_HEAT01.yaml") as fh:
        study = yaml.safe_load(fh)
    study["study_id"] = "BAD_VOCAB"
    study["comparisons"][0]["phenotype"] = "not_a_real_phenotype"
    path = tmp_path / "BAD_VOCAB.yaml"
    with open(path, "w") as fh:
        yaml.safe_dump(study, fh)

    result = validate_study_file(path)
    assert not result.valid
    assert any("phenotype" in e for e in result.errors)


def test_resilience_classification_mismatch_is_a_warning_not_an_error(tmp_path):
    with open("registry/studies/GIGAS_HEAT01.yaml") as fh:
        study = yaml.safe_load(fh)
    study["study_id"] = "WARN_MISMATCH"
    # survival's ontology default is "resilience"; declare something else to trigger a warning.
    study["comparisons"][0]["resilience_classification"] = "stress_response"
    path = tmp_path / "WARN_MISMATCH.yaml"
    with open(path, "w") as fh:
        yaml.safe_dump(study, fh)

    result = validate_study_file(path)
    assert result.valid
    assert result.warnings
