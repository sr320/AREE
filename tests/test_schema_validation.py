import glob
import textwrap

import yaml

from intake.schema_validate import validate_study_file


def test_all_demo_studies_are_valid():
    paths = sorted(glob.glob("registry/demo_studies/GIGAS_*.yaml"))
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
    with open("registry/demo_studies/GIGAS_HEAT01.yaml") as fh:
        study = yaml.safe_load(fh)
    study["study_id"] = "BAD_VOCAB"
    study["comparisons"][0]["phenotype"] = "not_a_real_phenotype"
    path = tmp_path / "BAD_VOCAB.yaml"
    with open(path, "w") as fh:
        yaml.safe_dump(study, fh)

    result = validate_study_file(path)
    assert not result.valid
    assert any("phenotype" in e for e in result.errors)


def test_batch_confounded_quality_flag_validates(tmp_path):
    with open("registry/demo_studies/GIGAS_HEAT01.yaml") as fh:
        study = yaml.safe_load(fh)
    study["study_id"] = "BATCH_FLAG"
    study["quality_flags"] = [*(study.get("quality_flags") or []), "batch_confounded"]
    path = tmp_path / "BATCH_FLAG.yaml"
    with open(path, "w") as fh:
        yaml.safe_dump(study, fh)

    result = validate_study_file(path)
    assert result.valid, result.errors


def test_resilience_classification_mismatch_is_a_warning_not_an_error(tmp_path):
    with open("registry/demo_studies/GIGAS_HEAT01.yaml") as fh:
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


def test_excluded_sample_missing_from_the_sample_sheet_is_rejected(tmp_path):
    """A typo in excluded_samples must not pass as an exclusion that never happens."""
    with open("registry/studies/CALLA2026_OSHV.yaml") as fh:
        study = yaml.safe_load(fh)
    study["comparisons"][0]["excluded_samples"][0]["sample_id"] = "DO_CT_44"
    path = tmp_path / "CALLA2026_OSHV.yaml"
    with open(path, "w") as fh:
        yaml.safe_dump(study, fh)

    result = validate_study_file(path)
    assert not result.valid
    assert any("DO_CT_44" in e for e in result.errors)


def test_excluded_sample_needs_a_reason_and_evidence(tmp_path):
    with open("registry/studies/CALLA2026_OSHV.yaml") as fh:
        study = yaml.safe_load(fh)
    del study["comparisons"][0]["excluded_samples"][0]["evidence"]
    path = tmp_path / "CALLA2026_OSHV.yaml"
    with open(path, "w") as fh:
        yaml.safe_dump(study, fh)

    assert not validate_study_file(path).valid


def _heat01_with(tmp_path, name, **comparison_fields):
    with open("registry/demo_studies/GIGAS_HEAT01.yaml") as fh:
        study = yaml.safe_load(fh)
    study["study_id"] = name
    study["comparisons"][0].update(comparison_fields)
    path = tmp_path / f"{name}.yaml"
    with open(path, "w") as fh:
        yaml.safe_dump(study, fh)
    return validate_study_file(path)


def test_meta_analysis_opt_out_with_a_reason_validates(tmp_path):
    reason = "Confounded field contrast: farm vs wild mixes disease with genetics and site."
    result = _heat01_with(tmp_path, "OPT_OUT", exclude_from_meta_analysis={"reason": reason})
    assert result.valid, result.errors


def test_meta_analysis_opt_out_needs_a_reason(tmp_path):
    assert not _heat01_with(tmp_path, "OPT_OUT_NO_REASON", exclude_from_meta_analysis={}).valid
    assert not _heat01_with(tmp_path, "OPT_OUT_BOOL", exclude_from_meta_analysis=True).valid


def test_meta_analysis_opt_out_cannot_also_be_primary(tmp_path):
    reason = "Confounded field contrast: farm vs wild mixes disease with genetics and site."
    both = _heat01_with(tmp_path, "OPT_OUT_PRIMARY", meta_analysis_primary=True,
                        exclude_from_meta_analysis={"reason": reason})
    assert not both.valid
    # meta_analysis_primary: false is just the default and does not conflict.
    assert _heat01_with(tmp_path, "OPT_OUT_NOT_PRIMARY", meta_analysis_primary=False,
                        exclude_from_meta_analysis={"reason": reason}).valid
