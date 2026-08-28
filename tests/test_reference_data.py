"""Tests for species identity and genome-assembly reference data.

These cover three things that were previously unenforced:

* `species` was free text, so an accepted synonym of the same animal would have
  split into two groups (and, once a second species is registered, two species
  could have pooled together);
* `data/reference/genome_assemblies.yaml` was read by nothing and shipped a
  placeholder accession;
* nothing recorded that identifiers are standardized against a *different*
  annotation from the assembly the source study used.
"""
from __future__ import annotations

import pytest
import yaml

from common import (
    GENOME_ASSEMBLIES_PATH,
    REPO_ROOT,
    load_assemblies,
    load_vocab,
    resolve_assembly,
    resolve_species,
)

# --------------------------------------------------------------------------- #
# Species identity
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "name",
    ["Magallana gigas", "Crassostrea gigas", "  crassostrea GIGAS  "],
)
def test_accepted_names_resolve_to_one_taxon(name):
    """Both genus names for the Pacific oyster must land on the same taxid."""
    resolved = resolve_species(name)
    assert resolved is not None
    assert resolved.ncbi_taxid == 29159
    assert resolved.scientific_name == "Magallana gigas"


def test_synonym_is_flagged_as_such():
    assert resolve_species("Crassostrea gigas").is_synonym is True
    assert resolve_species("Magallana gigas").is_synonym is False


def test_reported_name_is_preserved_not_overwritten():
    """AREE preserves what the source said; it does not silently rewrite it."""
    assert resolve_species("Crassostrea gigas").as_reported == "Crassostrea gigas"


@pytest.mark.parametrize("name", ["Ostrea edulis", "Homo sapiens", "", None, "gigas"])
def test_unknown_species_does_not_resolve(name):
    """A species must be added to the vocabulary deliberately, never guessed."""
    assert resolve_species(name) is None


def test_every_species_term_carries_a_taxid():
    for term in load_vocab("species").values():
        assert isinstance(term.get("ncbi_taxid"), int), term
        assert term.get("scientific_name")


# --------------------------------------------------------------------------- #
# Genome assemblies
# --------------------------------------------------------------------------- #


def test_no_placeholder_accessions_remain():
    """The reference file previously shipped PLACEHOLDER_CONFIRM_BEFORE_REAL_USE."""
    text = GENOME_ASSEMBLIES_PATH.read_text()
    assert "PLACEHOLDER" not in text.upper()


def test_assembly_accessions_are_well_formed():
    for record in load_assemblies().values():
        assert record["ncbi_assembly_accession"].startswith("GC"), record
        assert record.get("ncbi_taxid"), record


def test_exactly_one_assembly_is_the_ncbi_reference_per_species():
    by_species = {}
    for record in load_assemblies().values():
        if record.get("is_ncbi_reference"):
            by_species.setdefault(record["ncbi_taxid"], []).append(record["assembly_id"])
    for taxid, ids in by_species.items():
        assert len(ids) == 1, f"taxid {taxid} has multiple reference assemblies: {ids}"


@pytest.mark.parametrize(
    "value,expected",
    [
        ("cgigas_uk_roslin_v1", "cgigas_uk_roslin_v1"),
        ("GCA_902806645.1", "cgigas_uk_roslin_v1"),
        ("GCF_902806645.1", "cgigas_uk_roslin_v1"),
        # Studies commonly write the accession and the name together.
        ("GCA_902806645.1 (cgigas_uk_roslin_v1)", "cgigas_uk_roslin_v1"),
        ("GCF_963853765.1", "xbMagGiga1.1"),
        ("xbMagGiga1.1", "xbMagGiga1.1"),
    ],
)
def test_assembly_resolution(value, expected):
    assert resolve_assembly(value)["assembly_id"] == expected


@pytest.mark.parametrize("value", ["not_an_assembly", "", None, "GCA_000000000.1"])
def test_unknown_assembly_does_not_resolve(value):
    assert resolve_assembly(value) is None


# --------------------------------------------------------------------------- #
# Validation wiring
# --------------------------------------------------------------------------- #


def _study_fixture(tmp_path, **overrides):
    """A minimal valid study, written to a path whose stem matches its study_id."""
    from common import STUDIES_DIR

    base = yaml.safe_load((STUDIES_DIR / "HESSER2024_VCOR.yaml").read_text())
    base.update(overrides)
    path = tmp_path / f"{base['study_id']}.yaml"
    path.write_text(yaml.safe_dump(base))
    return path


def test_unknown_species_fails_validation(tmp_path):
    from intake.schema_validate import validate_study_file

    result = validate_study_file(
        _study_fixture(tmp_path, study_id="SPECIES_TEST", species="Ostrea edulis")
    )
    assert not result.valid
    assert any("species vocabulary" in e for e in result.errors)


def test_synonym_species_validates_with_a_warning(tmp_path):
    from intake.schema_validate import validate_study_file

    result = validate_study_file(
        _study_fixture(tmp_path, study_id="SYN_TEST", species="Crassostrea gigas")
    )
    assert result.valid
    assert any("accepted synonym" in w for w in result.warnings)


def test_canonical_species_validates_without_a_species_warning(tmp_path):
    from intake.schema_validate import validate_study_file

    result = validate_study_file(
        _study_fixture(tmp_path, study_id="CANON_TEST", species="Magallana gigas")
    )
    assert result.valid
    assert not any("synonym" in w for w in result.warnings)


def test_unknown_assembly_fails_validation(tmp_path):
    from intake.schema_validate import validate_study_file

    result = validate_study_file(
        _study_fixture(tmp_path, study_id="ASM_TEST", genome_assembly="my_local_assembly_v3")
    )
    assert not result.valid
    assert any("not a known assembly" in e for e in result.errors)


def test_all_registered_studies_pass_the_new_checks():
    """Every study already in the registry must satisfy the reference-data rules."""
    from common import STUDIES_DIR
    from intake.schema_validate import validate_study_file

    for path in sorted(STUDIES_DIR.glob("*.yaml")):
        if path.stem.startswith("_"):
            continue
        result = validate_study_file(path)
        assert result.valid, f"{path.name}: {result.errors}"


# --------------------------------------------------------------------------- #
# The annotation crossing must be visible in the evidence
# --------------------------------------------------------------------------- #


def test_real_crosswalk_reports_its_annotation_release():
    from harmonize.identifiers import crosswalk_annotation_release

    release = crosswalk_annotation_release(
        "data/reference/crosswalk/mgigas_gene_id_crosswalk.tsv"
    )
    assert release and "RS_2024_06" in release


def test_demo_crosswalk_reports_no_annotation_release():
    """The synthetic crosswalk has no provenance, and must not invent one."""
    from harmonize.identifiers import DEMO_CROSSWALK_PATH, crosswalk_annotation_release

    assert crosswalk_annotation_release(DEMO_CROSSWALK_PATH) is None


@pytest.fixture
def real_crosswalk():
    """Select the real crosswalk explicitly, as `make real-study` does via the env."""
    from harmonize.identifiers import set_crosswalk_path

    path = REPO_ROOT / "data/reference/crosswalk/mgigas_gene_id_crosswalk.tsv"
    if not path.exists():
        pytest.skip("real crosswalk not built; run `aree build-crosswalk`")
    set_crosswalk_path(path)
    yield path
    set_crosswalk_path(None)


def test_evidence_records_the_annotation_actually_used(real_crosswalk):
    """HESSER2024_VCOR declares the Roslin assembly but is standardized against
    the current NCBI annotation. That crossing must appear on every record."""
    from harmonize.core import harmonize_study

    df = harmonize_study("HESSER2024_VCOR", date_generated="2026-01-01")

    assert "cgigas_uk_roslin_v1" in df["genome_assembly"].iloc[0]
    releases = set(df["identifier_annotation_release"])
    assert len(releases) == 1
    assert "xbMagGiga1.1" in releases.pop()


def test_species_taxid_is_populated_for_every_record():
    from harmonize.core import harmonize_study

    df = harmonize_study("GIGAS_HEAT01", date_generated="2026-01-01")
    assert (df["species_taxid"] == 29159).all()
    # The reported name is preserved alongside the canonical id.
    assert (df["species"] == "Crassostrea gigas").all()


def test_meta_analysis_groups_on_canonical_species():
    """Grouping must use species_taxid, not the free-text name, so that the same
    animal under two accepted genus names is not split into two candidates."""
    from meta_analysis.run import GROUP_KEYS

    assert "species_taxid" in GROUP_KEYS
    assert "species" not in GROUP_KEYS


def test_harmonizing_a_study_still_marked_not_started_warns(real_crosswalk, tmp_path, monkeypatch):
    """analysis_status is curator-maintained and went stale once already."""
    import warnings as _warnings

    import harmonize.core as core

    study = core._load_study("HESSER2024_VCOR")
    study["analysis_status"] = "not_started"
    monkeypatch.setattr(core, "_load_study", lambda _sid: study)

    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        core.harmonize_study("HESSER2024_VCOR", date_generated="2026-01-01")

    assert any("analysis_status" in str(w.message) for w in caught)


def test_harmonizing_an_up_to_date_study_does_not_warn(real_crosswalk):
    import warnings as _warnings

    from harmonize.core import harmonize_study

    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        harmonize_study("HESSER2024_VCOR", date_generated="2026-01-01")

    assert not any("analysis_status" in str(w.message) for w in caught)
